"""
Translocation v3 — Same delta-delta interaction test for Intoxication and PA
=============================================================================
Identical method to make_translocation_v3_stats.py (AW version) but using
Intoxication and Protracted Abstinence replicates instead of AW.

For each condition:
  delta_CH[i] = CH_cond_rep[i] - mean(CH_naive)
  delta_SN[i] = SN_cond_rep[i] - mean(SN_naive)
  Welch t-test: delta_CH vs delta_SN
  BH FDR correction
  Positive interaction score = moved INTO Chromatin

Output: Translocation_IntoxPA_stats.xlsx
  Sheets per condition:
    All_proteins_Intox / All_proteins_PA
    Into_Chromatin_Intox / Into_Chromatin_PA
    G1_true_translocation_Intox / G1_true_translocation_PA   (CH up, SN down < -0.5)
    Into_SN_Intox / Into_SN_PA
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE       = os.path.join(SCRIPT_DIR, '../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_XLSX   = os.path.join(SCRIPT_DIR, 'Translocation_IntoxPA_stats.xlsx')

# ── Column groups ──────────────────────────────────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',   'Nuc_N-F-2',   'Nuc_N-F-3',   'Nuc_N-M-1',   'Nuc_N-M-2']

CH_I_COLS  = ['Chrom_I-F-1', 'Chrom_I-F-2', 'Chrom_I-F-3', 'Chrom_I-M-1', 'Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',   'Nuc_I-F-2',   'Nuc_I-F-3',   'Nuc_I-M-1',   'Nuc_I-M-2']

CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1',   'Nuc_AW-F-2',   'Nuc_AW-M-1',   'Nuc_AW-M-2']

CH_PA_COLS = ['Chrom_PA-F-1', 'Chrom_PA-F-2', 'Chrom_PA-F-3', 'Chrom_PA-M-1', 'Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1',   'Nuc_PA-F-2',   'Nuc_PA-F-3',   'Nuc_PA-M-1',   'Nuc_PA-M-2']

CONDITIONS = {
    'Intox': (CH_I_COLS,  SN_I_COLS),
    'PA':    (CH_PA_COLS, SN_PA_COLS),
}

# ── BH FDR ─────────────────────────────────────────────────────────────────────
def bh_correction(pvals):
    pvals = np.array(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    p_adj = pvals[order] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        p_adj[i] = min(p_adj[i], p_adj[i + 1])
    result = np.empty(n)
    result[order] = np.minimum(1.0, p_adj)
    return result

# ── Load sheets ────────────────────────────────────────────────────────────────
print('Loading Chromatin sheet (Keep + Review)...')
ch_raw = pd.read_excel(FILE, sheet_name='Chromatin')
ch_df  = ch_raw[ch_raw['Filter'].isin(['Keep', 'Review'])].copy().reset_index(drop=True)

print('Loading Soluble nuclear sheet (no filter)...')
sn_df  = pd.read_excel(FILE, sheet_name='Soluble nuclear').copy().reset_index(drop=True)

# Index by gene
for df in [ch_df, sn_df]:
    df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')

common_genes = sorted(set(ch_df.index) & set(sn_df.index))
print(f'Proteins in both fractions: {len(common_genes)}')

def get_vals(df, gene, cols):
    return pd.to_numeric(df.loc[gene, cols], errors='coerce').values \
           if gene in df.index else np.full(len(cols), np.nan)

def mean_fc(cond_cols, naive_cols, df, gene):
    c = get_vals(df, gene, cond_cols)
    n = get_vals(df, gene, naive_cols)
    return float(np.nanmean(c) - np.nanmean(n))

# ── Run analysis for each condition ───────────────────────────────────────────
all_results = {}

for cond_name, (ch_cond_cols, sn_cond_cols) in CONDITIONS.items():
    print(f'\n--- {cond_name} ---')
    rows = []
    for gene in common_genes:
        ch_naive = get_vals(ch_df, gene, CH_N_COLS)
        sn_naive = get_vals(sn_df, gene, SN_N_COLS)
        ch_cond  = get_vals(ch_df, gene, ch_cond_cols)
        sn_cond  = get_vals(sn_df, gene, sn_cond_cols)

        delta_ch = ch_cond - np.nanmean(ch_naive)
        delta_sn = sn_cond - np.nanmean(sn_naive)

        mean_dch = float(np.nanmean(delta_ch))
        mean_dsn = float(np.nanmean(delta_sn))
        interaction = mean_dch - mean_dsn

        d_ch = delta_ch[np.isfinite(delta_ch)]
        d_sn = delta_sn[np.isfinite(delta_sn)]
        if len(d_ch) >= 2 and len(d_sn) >= 2:
            _, p = stats.ttest_ind(d_ch, d_sn, equal_var=False)
        else:
            p = np.nan

        rows.append({
            'Gene':              gene,
            'Mean_delta_CH':     round(mean_dch, 4)      if np.isfinite(mean_dch)      else np.nan,
            'Mean_delta_SN':     round(mean_dsn, 4)      if np.isfinite(mean_dsn)      else np.nan,
            'Interaction_score': round(interaction, 4)   if np.isfinite(interaction)   else np.nan,
            # FC vs naive for all 3 conditions (for heatmap)
            'Chrom_FC_Intox': round(mean_fc(CH_I_COLS,  CH_N_COLS, ch_df, gene), 4),
            'Chrom_FC_AW':    round(mean_fc(CH_AW_COLS, CH_N_COLS, ch_df, gene), 4),
            'Chrom_FC_PA':    round(mean_fc(CH_PA_COLS, CH_N_COLS, ch_df, gene), 4),
            'SN_FC_Intox':    round(mean_fc(SN_I_COLS,  SN_N_COLS, sn_df, gene), 4),
            'SN_FC_AW':       round(mean_fc(SN_AW_COLS, SN_N_COLS, sn_df, gene), 4),
            'SN_FC_PA':       round(mean_fc(SN_PA_COLS, SN_N_COLS, sn_df, gene), 4),
            'p_value': p,
        })

    df_all = pd.DataFrame(rows)

    # BH FDR
    valid = df_all['p_value'].notna()
    df_all['p_adj_BH'] = np.nan
    df_all.loc[valid, 'p_adj_BH'] = bh_correction(df_all.loc[valid, 'p_value'].values)

    df_all['Direction']       = np.where(df_all['Interaction_score'] > 0, 'Into_Chromatin', 'Into_SN')
    df_all['Significant_p005'] = df_all['p_adj_BH'] < 0.05
    df_all = df_all.sort_values('Interaction_score', ascending=False).reset_index(drop=True)

    # Which FC column to use for SN down filter
    sn_fc_col = 'SN_FC_Intox' if cond_name == 'Intox' else 'SN_FC_PA'

    ch_sig  = df_all[(df_all['Direction'] == 'Into_Chromatin') & df_all['Significant_p005']].copy()
    sn_sig  = df_all[(df_all['Direction'] == 'Into_SN')        & df_all['Significant_p005']].copy()
    true_trans = ch_sig[ch_sig[sn_fc_col] < -0.5].copy()

    print(f'  Total tested:        {len(df_all)}')
    print(f'  Into Chromatin:      {len(ch_sig)}')
    print(f'  True translocation:  {len(true_trans)}  (CH up + SN down < -0.5)')
    print(f'  Into SN:             {len(sn_sig)}')

    all_results[cond_name] = {
        'all':   df_all,
        'ch':    ch_sig,
        'true':  true_trans,
        'sn':    sn_sig,
    }

# ── Save Excel ─────────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    for cond_name, res in all_results.items():
        res['all'].to_excel(writer,  sheet_name=f'All_proteins_{cond_name}',         index=False)
        res['ch'].to_excel(writer,   sheet_name=f'Into_Chromatin_{cond_name}',       index=False)
        res['true'].to_excel(writer, sheet_name=f'G1_true_translocation_{cond_name}',index=False)
        res['sn'].to_excel(writer,   sheet_name=f'Into_SN_{cond_name}',              index=False)

print(f'\nSaved: {OUT_XLSX}')
