"""
Translocation Union Heatmap — p < 0.10 threshold
==================================================
Same as make_translocation_union_heatmap.py but using BH FDR p_adj < 0.10
(~10% FDR, equivalent to −log2 ≥ 3.3 used in the per-fraction analyses).

Layout (top → bottom):
  • Significance annotation strip  (3 rows: Intox / AW / PA)
  • Chromatin FC band  (3 condition rows)
  • Soluble Nuclear FC band  (3 condition rows)

Excel output has Metascape-ready gene lists per condition × direction.

Outputs:
  Translocation_Union_stats_p10.xlsx
  Translocation_Union_heatmap_p10.pdf / .png
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE     = os.path.join(SCRIPT_DIR,
           '../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_Union_stats_p10.xlsx')
OUT_PDF  = os.path.join(SCRIPT_DIR, 'Translocation_Union_heatmap_p10.pdf')
OUT_PNG  = os.path.join(SCRIPT_DIR, 'Translocation_Union_heatmap_p10.png')

P_THRESH = 0.10   # ← relaxed from 0.05

# ── Column groups (AW-M-3 excluded) ──────────────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1','Chrom_N-F-2','Chrom_N-F-3','Chrom_N-M-1','Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',  'Nuc_N-F-2',  'Nuc_N-F-3',  'Nuc_N-M-1',  'Nuc_N-M-2']
CH_I_COLS  = ['Chrom_I-F-1','Chrom_I-F-2','Chrom_I-F-3','Chrom_I-M-1','Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',  'Nuc_I-F-2',  'Nuc_I-F-3',  'Nuc_I-M-1',  'Nuc_I-M-2']
CH_AW_COLS = ['Chrom_AW-F-1','Chrom_AW-F-2','Chrom_AW-M-1','Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1', 'Nuc_AW-F-2', 'Nuc_AW-M-1', 'Nuc_AW-M-2']
CH_PA_COLS = ['Chrom_PA-F-1','Chrom_PA-F-2','Chrom_PA-F-3','Chrom_PA-M-1','Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1', 'Nuc_PA-F-2', 'Nuc_PA-F-3', 'Nuc_PA-M-1', 'Nuc_PA-M-2']

CONDITIONS  = [('Intox', CH_I_COLS,  SN_I_COLS),
               ('AW',    CH_AW_COLS, SN_AW_COLS),
               ('PA',    CH_PA_COLS, SN_PA_COLS)]
COND_LABELS = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

def bh_correction(pvals):
    pvals = np.array(pvals, dtype=float)
    n = len(pvals); order = np.argsort(pvals)
    p_adj = pvals[order] * n / (np.arange(n) + 1)
    for i in range(n-2,-1,-1): p_adj[i] = min(p_adj[i], p_adj[i+1])
    result = np.empty(n); result[order] = np.minimum(1.0, p_adj)
    return result

# ── Load ──────────────────────────────────────────────────────────────────────
print('Loading data...')
ch_raw = pd.read_excel(FILE, sheet_name='Chromatin')
ch_df  = ch_raw[ch_raw['Filter'].isin(['Keep','Review'])].copy().reset_index(drop=True)
sn_df  = pd.read_excel(FILE, sheet_name='Soluble nuclear').copy().reset_index(drop=True)
for df in [ch_df, sn_df]:
    df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
common_genes = sorted(set(ch_df.index) & set(sn_df.index))
print(f'Proteins in both fractions: {len(common_genes)}')

def get_vals(df, gene, cols):
    return pd.to_numeric(df.loc[gene, cols], errors='coerce').values \
           if gene in df.index else np.full(len(cols), np.nan)

def mean_fc_val(df, gene, cond_cols, naive_cols):
    c = get_vals(df, gene, cond_cols); n = get_vals(df, gene, naive_cols)
    return float(np.nanmean(c) - np.nanmean(n))

# ── Delta-delta test ──────────────────────────────────────────────────────────
print('Running unified interaction test...')
rows = []
for gene in common_genes:
    ch_naive = get_vals(ch_df, gene, CH_N_COLS)
    sn_naive = get_vals(sn_df, gene, SN_N_COLS)
    row = {'Gene': gene}
    for cond_name, ch_cols, sn_cols in CONDITIONS:
        ch_cond  = get_vals(ch_df, gene, ch_cols)
        sn_cond  = get_vals(sn_df, gene, sn_cols)
        delta_ch = ch_cond - np.nanmean(ch_naive)
        delta_sn = sn_cond - np.nanmean(sn_naive)
        interaction = float(np.nanmean(delta_ch) - np.nanmean(delta_sn))
        d_ch = delta_ch[np.isfinite(delta_ch)]; d_sn = delta_sn[np.isfinite(delta_sn)]
        p = stats.ttest_ind(d_ch, d_sn, equal_var=False)[1] \
            if len(d_ch) >= 2 and len(d_sn) >= 2 else np.nan
        row[f'Interaction_{cond_name}'] = round(interaction, 4) if np.isfinite(interaction) else np.nan
        row[f'p_{cond_name}']           = p
        row[f'Chrom_FC_{cond_name}']    = round(mean_fc_val(ch_df, gene, ch_cols, CH_N_COLS), 4)
        row[f'SN_FC_{cond_name}']       = round(mean_fc_val(sn_df, gene, sn_cols, SN_N_COLS), 4)
    rows.append(row)

df_all = pd.DataFrame(rows)

# BH FDR per condition
for cond_name, _, _ in CONDITIONS:
    valid = df_all[f'p_{cond_name}'].notna()
    df_all[f'p_adj_{cond_name}'] = np.nan
    df_all.loc[valid, f'p_adj_{cond_name}'] = bh_correction(
        df_all.loc[valid, f'p_{cond_name}'].values)
    df_all[f'Sig_{cond_name}'] = df_all[f'p_adj_{cond_name}'] < P_THRESH
    df_all[f'Direction_{cond_name}'] = np.where(
        df_all[f'Sig_{cond_name}'] & (df_all[f'Interaction_{cond_name}'] > 0), 'Into_Chromatin',
        np.where(
        df_all[f'Sig_{cond_name}'] & (df_all[f'Interaction_{cond_name}'] < 0), 'Into_SN',
        'NS'))

print('\n── Per-condition summary ──')
for cond_name, _, _ in CONDITIONS:
    n_ch = (df_all[f'Direction_{cond_name}'] == 'Into_Chromatin').sum()
    n_sn = (df_all[f'Direction_{cond_name}'] == 'Into_SN').sum()
    print(f'  {cond_name}: Into_Chromatin={n_ch}  Into_SN={n_sn}')

union_mask = df_all['Sig_Intox'] | df_all['Sig_AW'] | df_all['Sig_PA']
df_union   = df_all[union_mask].copy().reset_index(drop=True)
print(f'\nUnion proteins (sig in ≥1 condition): {len(df_union)}')

# ── Save Excel (Metascape-ready) ──────────────────────────────────────────────
# Build flat gene-list sheet for easy Metascape copy-paste
meta_rows = []
for cond_name, _, _ in CONDITIONS:
    for dirn in ['Into_Chromatin', 'Into_SN']:
        mask = df_all[f'Direction_{cond_name}'] == dirn
        genes = df_all.loc[mask, 'Gene'].tolist()
        label = f'{cond_name}_{dirn}'
        for g in genes:
            meta_rows.append({'List': label, 'Condition': cond_name,
                               'Direction': dirn, 'Gene': g})
df_meta = pd.DataFrame(meta_rows)

with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    df_all.to_excel(writer,   sheet_name='All_proteins',   index=False)
    df_union.to_excel(writer, sheet_name='Union_sig',      index=False)
    df_meta.to_excel(writer,  sheet_name='Metascape_lists', index=False)
    for cond_name, _, _ in CONDITIONS:
        for dirn, short in [('Into_Chromatin','Into_CH'), ('Into_SN','Into_SN')]:
            sub = df_all[df_all[f'Direction_{cond_name}'] == dirn][
                ['Gene', f'Chrom_FC_{cond_name}', f'SN_FC_{cond_name}',
                 f'Interaction_{cond_name}', f'p_adj_{cond_name}',
                 f'Direction_{cond_name}']].copy()
            sub.columns = ['Gene','Chrom_FC','SN_FC','Interaction','p_adj','Direction']
            sub.to_excel(writer, sheet_name=f'{short}_{cond_name}', index=False)
print(f'Saved: {OUT_XLSX}')

# ── Sort: Into_Chromatin (positive sum) → Into_SN (negative sum) ──────────────
ti_cols = [f'Interaction_{c}' for c, _, _ in CONDITIONS]
df_union['sum_interaction'] = df_union[ti_cols].fillna(0).sum(axis=1)
df_plot = df_union.sort_values('sum_interaction', ascending=False).reset_index(drop=True)
n_prot  = len(df_plot)
print(f'Sorted by sum interaction. n = {n_prot}')

# ── Colormaps ─────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0','#2B7FD4','#6AAEE0','#AECFE8','#DDEEF8',
     '#FFFFFF',
     '#FDD8E7','#F5A0BC','#EE5F8B','#E8305A','#C01E42'], N=512)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)
ANNOT_CMAP = ListedColormap(['#2B7FD4', '#E8E8E8', '#E8305A'])
ANNOT_NORM = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ANNOT_CMAP.N)

# ── Layout ────────────────────────────────────────────────────────────────────
INCHES_PER_PROT = 0.007
ROW_H    = 0.50; N_COND = 3; BAND_H = N_COND * ROW_H
ANNOT_ROW_H = 0.12; ANNOT_H = N_COND * ANNOT_ROW_H
LEFT_IN = 2.50; RIGHT_IN = 2.80; TOP_PAD = 0.55; BOT_PAD = 1.20
ANNOT_GAP = 0.14; SECTION_GAP = 0.32
CBAR_W_IN = 2.8; CBAR_H_IN = 0.14; CBAR_BOT = 0.48
COMP_BG = {'Chromatin': '#FBF9F0', 'Soluble Nuclear': '#F0F6FB'}

hm_w  = max(n_prot * INCHES_PER_PROT, 3.0)
FIG_W = LEFT_IN + hm_w + RIGHT_IN
FIG_H = TOP_PAD + ANNOT_H + ANNOT_GAP + BAND_H + SECTION_GAP + BAND_H + BOT_PAD

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
y_cursor = FIG_H - TOP_PAD

# Annotation strip
annot_mat = np.zeros((3, n_prot), dtype=float)
for ci, (cond_name, _, _) in enumerate(CONDITIONS):
    for pi, (_, row) in enumerate(df_plot.iterrows()):
        d = row[f'Direction_{cond_name}']
        annot_mat[ci, pi] = 1 if d == 'Into_Chromatin' else (-1 if d == 'Into_SN' else 0)

y_annot_bot = y_cursor - ANNOT_H
ax_annot = fig.add_axes([LEFT_IN/FIG_W, y_annot_bot/FIG_H, hm_w/FIG_W, ANNOT_H/FIG_H])
ax_annot.imshow(annot_mat, aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM, interpolation='nearest')
for r in range(1, N_COND): ax_annot.axhline(r-0.5, color='white', lw=0.8, zorder=3)
ax_annot.set_xticks([]); ax_annot.set_yticks([])
for sp in ax_annot.spines.values(): sp.set_visible(False)
for r, cname in enumerate(COND_LABELS):
    y_frac = (y_annot_bot + ANNOT_H - (r+0.5)*ANNOT_ROW_H) / FIG_H
    fig.text((LEFT_IN+hm_w+0.12)/FIG_W, y_frac, cname, ha='left', va='center',
             fontsize=6.5, fontfamily='Arial', color='#333333')
fig.text((LEFT_IN-0.12)/FIG_W, (y_annot_bot+ANNOT_H/2)/FIG_H,
         'Significance\n& Direction', ha='right', va='center',
         fontsize=7, fontfamily='Arial', color='#555555', style='italic')
y_cursor = y_annot_bot - ANNOT_GAP

# Two-band heatmap
for label, fc_cols in [
    ('Chromatin',       ['Chrom_FC_Intox', 'Chrom_FC_AW', 'Chrom_FC_PA']),
    ('Soluble Nuclear', ['SN_FC_Intox',    'SN_FC_AW',    'SN_FC_PA']),
]:
    y_band_bot = y_cursor - BAND_H
    mat = np.vstack([df_plot[col].fillna(0).values for col in fc_cols])
    ax  = fig.add_axes([LEFT_IN/FIG_W, y_band_bot/FIG_H, hm_w/FIG_W, BAND_H/FIG_H],
                       facecolor=COMP_BG[label])
    im  = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    for r in range(1, N_COND): ax.axhline(r-0.5, color='white', lw=1.2, zorder=3)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    for r, cname in enumerate(COND_LABELS):
        y_frac = (y_band_bot + BAND_H - (r+0.5)*ROW_H) / FIG_H
        fig.text((LEFT_IN+hm_w+0.12)/FIG_W, y_frac, cname, ha='left', va='center',
                 fontsize=8, fontfamily='Arial', color='#222222')
    mid_y = (y_band_bot + BAND_H/2) / FIG_H
    fig.text((LEFT_IN-0.12)/FIG_W, mid_y, label, ha='right', va='center',
             fontsize=10, fontweight='bold', fontfamily='Arial', color='#111111')
    fig.text((LEFT_IN-0.12)/FIG_W, mid_y-0.055, f'n = {n_prot}', ha='right', va='top',
             fontsize=7.5, fontfamily='Arial', color='#555555', style='italic')
    y_cursor = y_band_bot - SECTION_GAP

# Legend
legend_patches = [
    mpatches.Patch(color='#E8305A', label='Into Chromatin'),
    mpatches.Patch(color='#2B7FD4', label='Into SN'),
    mpatches.Patch(color='#E8E8E8', label='NS', linewidth=0.5, edgecolor='#AAAAAA'),
]
fig.legend(handles=legend_patches, loc='lower right',
           bbox_to_anchor=(0.98, CBAR_BOT/FIG_H+0.04),
           fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
           title='Translocation', title_fontsize=7)

# Colorbar
cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
ax_cb = fig.add_axes([cbar_left, CBAR_BOT/FIG_H, CBAR_W_IN/FIG_W, CBAR_H_IN/FIG_H])
cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
cb.set_ticks([-3,-2,-1,0,1,2,3]); cb.ax.tick_params(labelsize=7)
ax_cb.set_xlabel('Log2 Fold Change vs Naive', fontsize=8, fontfamily='Arial',
                 color='#444444', labelpad=4)

# Title & footnote
fig.text(0.5, (FIG_H-TOP_PAD*0.35)/FIG_H,
         f'Protein Translocation — Union Across All Conditions  (n = {n_prot})  [p < 0.10]',
         ha='center', va='top', fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')
fig.text(0.5, 0.008,
         'Union: BH FDR p_adj < 0.10 (~10% FDR) in ≥1 condition  |  '
         'Hierarchical clustering (Ward) on translocation index  |  AW-M-3 excluded  |  Chromatin: Keep+Review only',
         ha='center', va='bottom', fontsize=5.5, fontfamily='Arial', color='#888888')

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
