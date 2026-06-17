"""
Translocation v3 — Proteins that moved INTO Chromatin during AW (vs Naive)
===========================================================================
Test (per protein):
  delta_CH  = CH_AW_rep  - mean(CH_Naive)   [4 values]
  delta_SN  = SN_AW_rep  - mean(SN_Naive)   [4 values]
  Welch's t-test: delta_CH vs delta_SN
  → significant + delta_CH > delta_SN  means the protein
    increased MORE in chromatin than in SN during AW

Settings
--------
  Chromatin : Keep + Review filter, AW-M-3 excluded (absent from file)
  SN        : no filter, AW-M-3 excluded (absent from file)
  Source    : EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx

Output
------
  Translocation_v3_stats.xlsx   (3 sheets: All_proteins, CH_higher_sig, SN_higher_sig)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE       = os.path.join(SCRIPT_DIR, '../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_XLSX   = os.path.join(SCRIPT_DIR, 'Translocation_v3_stats.xlsx')

# ── Load sheets ────────────────────────────────────────────────────────────────
print('Loading Chromatin sheet (Keep + Review)…')
ch_raw = pd.read_excel(FILE, sheet_name='Chromatin')
ch_df  = ch_raw[ch_raw['Filter'].isin(['Keep', 'Review'])].copy().reset_index(drop=True)

print('Loading Soluble nuclear sheet (no filter)…')
sn_df  = pd.read_excel(FILE, sheet_name='Soluble nuclear').copy().reset_index(drop=True)

# ── Rep columns ────────────────────────────────────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',   'Nuc_N-F-2',   'Nuc_N-F-3',   'Nuc_N-M-1',   'Nuc_N-M-2']
CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1',   'Nuc_AW-F-2',   'Nuc_AW-M-1',   'Nuc_AW-M-2']

# Also load Intox and PA cols for heatmap FC values
CH_I_COLS  = ['Chrom_I-F-1', 'Chrom_I-F-2', 'Chrom_I-F-3', 'Chrom_I-M-1', 'Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',   'Nuc_I-F-2',   'Nuc_I-F-3',   'Nuc_I-M-1',   'Nuc_I-M-2']
CH_PA_COLS = ['Chrom_PA-F-1', 'Chrom_PA-F-2', 'Chrom_PA-F-3', 'Chrom_PA-M-1', 'Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1',   'Nuc_PA-F-2',   'Nuc_PA-F-3',   'Nuc_PA-M-1',   'Nuc_PA-M-2']

# ── BH FDR correction ──────────────────────────────────────────────────────────
def bh_correction(pvals):
    pvals = np.array(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    p_adj_sorted = pvals[order] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        p_adj_sorted[i] = min(p_adj_sorted[i], p_adj_sorted[i + 1])
    result = np.empty(n)
    result[order] = np.minimum(1.0, p_adj_sorted)
    return result

# ── Index by Gene symbol ───────────────────────────────────────────────────────
ch_df['Gene symbol'] = ch_df['Gene symbol'].astype(str).str.strip()
sn_df['Gene symbol'] = sn_df['Gene symbol'].astype(str).str.strip()
ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')

common_genes = sorted(set(ch_df.index) & set(sn_df.index))
print(f'Proteins in both fractions: {len(common_genes)}')

def get_vals(df, gene, cols):
    return pd.to_numeric(df.loc[gene, cols], errors='coerce').values if gene in df.index else np.full(len(cols), np.nan)

def mean_fc(cond_cols, naive_cols, df, gene):
    c = get_vals(df, gene, cond_cols)
    n = get_vals(df, gene, naive_cols)
    return float(np.nanmean(c) - np.nanmean(n))

# ── Test per protein ───────────────────────────────────────────────────────────
rows = []
for gene in common_genes:
    ch_naive = get_vals(ch_df, gene, CH_N_COLS)
    sn_naive = get_vals(sn_df, gene, SN_N_COLS)
    ch_aw    = get_vals(ch_df, gene, CH_AW_COLS)
    sn_aw    = get_vals(sn_df, gene, SN_AW_COLS)

    # delta = AW rep − mean(naive)  → how much did each AW animal deviate from naive
    delta_ch = ch_aw - np.nanmean(ch_naive)
    delta_sn = sn_aw - np.nanmean(sn_naive)

    mean_delta_ch = float(np.nanmean(delta_ch))
    mean_delta_sn = float(np.nanmean(delta_sn))
    interaction   = mean_delta_ch - mean_delta_sn   # positive = moved into chromatin

    d_ch = delta_ch[np.isfinite(delta_ch)]
    d_sn = delta_sn[np.isfinite(delta_sn)]
    if len(d_ch) >= 2 and len(d_sn) >= 2:
        _, p = stats.ttest_ind(d_ch, d_sn, equal_var=False)
    else:
        p = np.nan

    rows.append({
        'Gene':              gene,
        'Filter_Chrom':      ch_raw.loc[ch_raw['Gene symbol'].astype(str).str.strip() == gene, 'Filter'].values[0]
                             if gene in ch_raw['Gene symbol'].astype(str).str.strip().values else np.nan,
        # Raw AW reps
        'Chrom_AW_F1': ch_df.loc[gene, 'Chrom_AW-F-1'] if gene in ch_df.index else np.nan,
        'Chrom_AW_F2': ch_df.loc[gene, 'Chrom_AW-F-2'] if gene in ch_df.index else np.nan,
        'Chrom_AW_M1': ch_df.loc[gene, 'Chrom_AW-M-1'] if gene in ch_df.index else np.nan,
        'Chrom_AW_M2': ch_df.loc[gene, 'Chrom_AW-M-2'] if gene in ch_df.index else np.nan,
        'SN_AW_F1': sn_df.loc[gene, 'Nuc_AW-F-1'] if gene in sn_df.index else np.nan,
        'SN_AW_F2': sn_df.loc[gene, 'Nuc_AW-F-2'] if gene in sn_df.index else np.nan,
        'SN_AW_M1': sn_df.loc[gene, 'Nuc_AW-M-1'] if gene in sn_df.index else np.nan,
        'SN_AW_M2': sn_df.loc[gene, 'Nuc_AW-M-2'] if gene in sn_df.index else np.nan,
        # Delta from naive (per fraction)
        'Mean_delta_CH': round(mean_delta_ch, 4) if np.isfinite(mean_delta_ch) else np.nan,
        'Mean_delta_SN': round(mean_delta_sn, 4) if np.isfinite(mean_delta_sn) else np.nan,
        'Interaction_score': round(interaction, 4) if np.isfinite(interaction) else np.nan,
        # Log2 FC vs naive per condition (for heatmap)
        'Chrom_FC_Intox': round(mean_fc(CH_I_COLS,  CH_N_COLS, ch_df, gene), 4),
        'Chrom_FC_AW':    round(mean_fc(CH_AW_COLS, CH_N_COLS, ch_df, gene), 4),
        'Chrom_FC_PA':    round(mean_fc(CH_PA_COLS, CH_N_COLS, ch_df, gene), 4),
        'SN_FC_Intox':    round(mean_fc(SN_I_COLS,  SN_N_COLS, sn_df, gene), 4),
        'SN_FC_AW':       round(mean_fc(SN_AW_COLS, SN_N_COLS, sn_df, gene), 4),
        'SN_FC_PA':       round(mean_fc(SN_PA_COLS, SN_N_COLS, sn_df, gene), 4),
        'p_value': p,
    })

df_all = pd.DataFrame(rows)

# ── BH FDR ─────────────────────────────────────────────────────────────────────
valid = df_all['p_value'].notna()
df_all['p_adj_BH'] = np.nan
df_all.loc[valid, 'p_adj_BH'] = bh_correction(df_all.loc[valid, 'p_value'].values)

df_all['Direction'] = np.where(df_all['Interaction_score'] > 0, 'Into_Chromatin', 'Into_SN')
df_all['Significant_p005'] = df_all['p_adj_BH'] < 0.05
df_all = df_all.sort_values('Interaction_score', ascending=False).reset_index(drop=True)

ch_higher_sig = df_all[(df_all['Direction'] == 'Into_Chromatin') & df_all['Significant_p005']].copy()
sn_higher_sig = df_all[(df_all['Direction'] == 'Into_SN')        & df_all['Significant_p005']].copy()

print(f'\nResults:')
print(f'  Total proteins tested        : {len(df_all)}')
print(f'  Into Chromatin (p_adj<0.05)  : {len(ch_higher_sig)}')
print(f'  Into SN        (p_adj<0.05)  : {len(sn_higher_sig)}')

with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    df_all.to_excel(writer,        sheet_name='All_proteins',    index=False)
    ch_higher_sig.to_excel(writer, sheet_name='Into_Chromatin',  index=False)
    sn_higher_sig.to_excel(writer, sheet_name='Into_SN',         index=False)

print(f'\nSaved: {OUT_XLSX}')
