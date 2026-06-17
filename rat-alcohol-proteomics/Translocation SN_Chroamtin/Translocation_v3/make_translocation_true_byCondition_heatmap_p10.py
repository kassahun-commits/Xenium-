"""
Translocation TRUE Union Heatmap — by-condition layout, p < 0.10
================================================================
True translocation proteins (p_adj < 0.10), shown as 3 condition bands
with Chromatin and SN as adjacent rows per band.

True translocation criteria:
  Into Chromatin: Sig + Into_Chromatin + SN_FC < -0.5
  Into SN:        Sig + Into_SN + SN_FC > 0 + Chrom_FC < -0.5

Outputs:
  Translocation_TrueUnion_byCondition_heatmap_p10.pdf / .png
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
OUT_PDF  = os.path.join(SCRIPT_DIR, 'Translocation_TrueUnion_byCondition_heatmap_p10.pdf')
OUT_PNG  = os.path.join(SCRIPT_DIR, 'Translocation_TrueUnion_byCondition_heatmap_p10.png')

P_THRESH  = 0.10
FC_THRESH = 0.5

# ── Column groups ─────────────────────────────────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1','Chrom_N-F-2','Chrom_N-F-3','Chrom_N-M-1','Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',  'Nuc_N-F-2',  'Nuc_N-F-3',  'Nuc_N-M-1',  'Nuc_N-M-2']
CH_I_COLS  = ['Chrom_I-F-1','Chrom_I-F-2','Chrom_I-F-3','Chrom_I-M-1','Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',  'Nuc_I-F-2',  'Nuc_I-F-3',  'Nuc_I-M-1',  'Nuc_I-M-2']
CH_AW_COLS = ['Chrom_AW-F-1','Chrom_AW-F-2','Chrom_AW-M-1','Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1', 'Nuc_AW-F-2', 'Nuc_AW-M-1', 'Nuc_AW-M-2']
CH_PA_COLS = ['Chrom_PA-F-1','Chrom_PA-F-2','Chrom_PA-F-3','Chrom_PA-M-1','Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1', 'Nuc_PA-F-2', 'Nuc_PA-F-3', 'Nuc_PA-M-1', 'Nuc_PA-M-2']

CONDITIONS  = [('Intox', CH_I_COLS, SN_I_COLS),
               ('AW',    CH_AW_COLS, SN_AW_COLS),
               ('PA',    CH_PA_COLS, SN_PA_COLS)]
COND_LABELS = ['Intoxication', 'Acute\nWithdrawal', 'Protracted\nAbstinence']

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

def get_vals(df, gene, cols):
    return pd.to_numeric(df.loc[gene, cols], errors='coerce').values \
           if gene in df.index else np.full(len(cols), np.nan)

def mean_fc(df, gene, cond_cols, naive_cols):
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
        row[f'Chrom_FC_{cond_name}']    = round(mean_fc(ch_df, gene, ch_cols, CH_N_COLS), 4)
        row[f'SN_FC_{cond_name}']       = round(mean_fc(sn_df, gene, sn_cols, SN_N_COLS), 4)
    rows.append(row)

df_all = pd.DataFrame(rows)

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

# True translocation filter
true_ch_any = pd.Series(False, index=df_all.index)
true_sn_any = pd.Series(False, index=df_all.index)
for cond_name, _, _ in CONDITIONS:
    sig   = df_all[f'Sig_{cond_name}']
    dirn  = df_all[f'Direction_{cond_name}']
    sn_fc = pd.to_numeric(df_all[f'SN_FC_{cond_name}'],    errors='coerce')
    ch_fc = pd.to_numeric(df_all[f'Chrom_FC_{cond_name}'], errors='coerce')
    true_ch = sig & (dirn == 'Into_Chromatin') & (sn_fc < -FC_THRESH)
    true_sn = sig & (dirn == 'Into_SN')        & (sn_fc >  0)          & (ch_fc < -FC_THRESH)
    true_ch_any = true_ch_any | true_ch
    true_sn_any = true_sn_any | true_sn
    print(f'{cond_name}: True Into_Chromatin={true_ch.sum()}  True Into_SN={true_sn.sum()}')

df_true = df_all[true_ch_any | true_sn_any].copy().reset_index(drop=True)
print(f'\nTrue union: {len(df_true)} proteins')

# ── Sort: Into_Chromatin (positive sum) → Into_SN (negative sum) ──────────────
ti_cols = [f'Interaction_{c}' for c, _, _ in CONDITIONS]
df_true['sum_interaction'] = df_true[ti_cols].fillna(0).sum(axis=1)
df_plot = df_true.sort_values('sum_interaction', ascending=False).reset_index(drop=True)
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
FRAC_ROW_H  = 0.50
BAND_H      = 2 * FRAC_ROW_H
ANNOT_ROW_H = 0.13; ANNOT_H = 3 * ANNOT_ROW_H
ANNOT_GAP   = 0.18; SECTION_GAP = 0.28; N_COND = 3
LEFT_IN = 2.60; RIGHT_IN = 2.80; TOP_PAD = 0.55; BOT_PAD = 1.20
CBAR_W_IN = 2.8; CBAR_H_IN = 0.14; CBAR_BOT = 0.46

hm_w  = max(n_prot * 0.007, 3.0)
FIG_W = LEFT_IN + hm_w + RIGHT_IN
FIG_H = (TOP_PAD + ANNOT_H + ANNOT_GAP
         + N_COND * BAND_H + (N_COND-1) * SECTION_GAP + BOT_PAD)

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
    fig.text((LEFT_IN+hm_w+0.10)/FIG_W, y_frac, cname, ha='left', va='center',
             fontsize=6, fontfamily='Arial', color='#333333')
fig.text((LEFT_IN-0.12)/FIG_W, (y_annot_bot+ANNOT_H/2)/FIG_H,
         'Significance\n& Direction', ha='right', va='center',
         fontsize=7, fontfamily='Arial', color='#555555', style='italic')
y_cursor = y_annot_bot - ANNOT_GAP

# Three condition bands
COMP_BG = {'Chromatin': '#FBF9F0', 'SN': '#F0F6FB'}

for ci, (cond_name, _, _) in enumerate(CONDITIONS):
    y_band_bot = y_cursor - BAND_H

    ch_mat = df_plot[f'Chrom_FC_{cond_name}'].fillna(0).values.reshape(1, -1)
    ax_ch  = fig.add_axes([LEFT_IN/FIG_W, (y_band_bot+FRAC_ROW_H)/FIG_H,
                            hm_w/FIG_W, FRAC_ROW_H/FIG_H],
                           facecolor=COMP_BG['Chromatin'])
    im = ax_ch.imshow(ch_mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_ch.set_xticks([]); ax_ch.set_yticks([])
    for sp in ax_ch.spines.values(): sp.set_visible(False)
    fig.text((LEFT_IN+hm_w+0.10)/FIG_W,
             (y_band_bot+FRAC_ROW_H+FRAC_ROW_H/2)/FIG_H,
             'Chromatin', ha='left', va='center',
             fontsize=7.5, fontfamily='Arial', color='#111111')

    sn_mat = df_plot[f'SN_FC_{cond_name}'].fillna(0).values.reshape(1, -1)
    ax_sn  = fig.add_axes([LEFT_IN/FIG_W, y_band_bot/FIG_H,
                            hm_w/FIG_W, FRAC_ROW_H/FIG_H],
                           facecolor=COMP_BG['SN'])
    ax_sn.imshow(sn_mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_sn.set_xticks([]); ax_sn.set_yticks([])
    for sp in ax_sn.spines.values(): sp.set_visible(False)
    fig.text((LEFT_IN+hm_w+0.10)/FIG_W, (y_band_bot+FRAC_ROW_H/2)/FIG_H,
             'Soluble\nNuclear', ha='left', va='center',
             fontsize=7.5, fontfamily='Arial', color='#111111')

    mid_y = (y_band_bot + BAND_H/2) / FIG_H
    fig.text((LEFT_IN-0.12)/FIG_W, mid_y,
             COND_LABELS[ci].replace('\n', ' '),
             ha='right', va='center',
             fontsize=10, fontweight='bold', fontfamily='Arial', color='#111111')

    sep_y = (y_band_bot + FRAC_ROW_H) / FIG_H
    fig.add_artist(plt.Line2D(
        [LEFT_IN/FIG_W, (LEFT_IN+hm_w)/FIG_W], [sep_y, sep_y],
        transform=fig.transFigure, color='#CCCCCC', linewidth=0.8, zorder=5))

    y_cursor = y_band_bot - SECTION_GAP

# Legend
legend_patches = [
    mpatches.Patch(color='#E8305A', label='Into Chromatin'),
    mpatches.Patch(color='#2B7FD4', label='Into SN'),
    mpatches.Patch(color='#E8E8E8', label='NS', linewidth=0.5, edgecolor='#AAAAAA'),
]
fig.legend(handles=legend_patches, loc='lower right',
           bbox_to_anchor=(0.98, CBAR_BOT/FIG_H+0.045),
           fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
           title='Translocation\n(annotation strip)', title_fontsize=7)

# Colorbar
cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
ax_cb = fig.add_axes([cbar_left, CBAR_BOT/FIG_H, CBAR_W_IN/FIG_W, CBAR_H_IN/FIG_H])
cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
cb.set_ticks([-3,-2,-1,0,1,2,3]); cb.ax.tick_params(labelsize=7)
ax_cb.set_xlabel('Log2 Fold Change vs Naive', fontsize=8, fontfamily='Arial',
                 color='#444444', labelpad=4)

# Title & footnote
fig.text(0.5, (FIG_H-TOP_PAD*0.35)/FIG_H,
         f'True Protein Translocation — Chromatin & SN by Condition  (n = {n_prot})  [p < 0.10]',
         ha='center', va='top', fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')
fig.text(0.5, 0.008,
         f'True Into_CH: Sig(p<0.10) + Into_Chromatin + SN_FC<-0.5  |  '
         f'True Into_SN: Sig(p<0.10) + Into_SN + SN_FC>0 + Chrom_FC<-0.5  |  '
         f'Union ≥1 condition  |  Ward clustering  |  AW-M-3 excluded',
         ha='center', va='bottom', fontsize=5.5, fontfamily='Arial', color='#888888')

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
