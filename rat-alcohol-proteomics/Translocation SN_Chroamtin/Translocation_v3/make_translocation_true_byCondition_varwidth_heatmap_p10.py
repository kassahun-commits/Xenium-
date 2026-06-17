"""
Translocation TRUE — by-condition, variable-width heatmap, p < 0.10
====================================================================
Each condition band shows ONLY the proteins with true translocation in
THAT condition, sorted by interaction value (Into_Chromatin → Into_SN).
Band width is proportional to protein count, making the n differences
immediately visible (AW >> Intox >> PA).

All bands are left-aligned; the narrower Intox and PA bands end before
the right margin, leaving white space that honestly represents scale.

Outputs:
  Translocation_TrueUnion_byCondition_varwidth_heatmap_p10.pdf / .png
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
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
OUT_PDF  = os.path.join(SCRIPT_DIR, 'Translocation_TrueUnion_byCondition_varwidth_heatmap_p10.pdf')
OUT_PNG  = os.path.join(SCRIPT_DIR, 'Translocation_TrueUnion_byCondition_varwidth_heatmap_p10.png')

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

# ── Per-condition true translocation protein lists ────────────────────────────
IPP       = 0.012   # inches per protein column
MIN_COL_W = 0.40    # minimum band width (inches) so PA is visible

cond_data = {}
for cond_name, _, _ in CONDITIONS:
    sn_fc = pd.to_numeric(df_all[f'SN_FC_{cond_name}'],    errors='coerce')
    ch_fc = pd.to_numeric(df_all[f'Chrom_FC_{cond_name}'], errors='coerce')
    sig   = df_all[f'Sig_{cond_name}']
    dirn  = df_all[f'Direction_{cond_name}']

    true_ch = sig & (dirn == 'Into_Chromatin') & (sn_fc < -FC_THRESH)
    true_sn = sig & (dirn == 'Into_SN')        & (sn_fc >  0)    & (ch_fc < -FC_THRESH)

    df_cond = df_all[true_ch | true_sn].copy()
    # True direction label per protein for this condition
    df_cond['true_dirn'] = np.where(
        true_ch.loc[df_cond.index], 'Into_Chromatin', 'Into_SN')
    # Sort: most positive interaction first (Into_CH), most negative last (Into_SN)
    df_cond = df_cond.sort_values(
        f'Interaction_{cond_name}', ascending=False).reset_index(drop=True)

    n = len(df_cond)
    w = max(n * IPP, MIN_COL_W)
    n_ch = (df_cond['true_dirn'] == 'Into_Chromatin').sum()
    n_sn = (df_cond['true_dirn'] == 'Into_SN').sum()
    cond_data[cond_name] = {'df': df_cond, 'n': n, 'w': w,
                             'n_ch': n_ch, 'n_sn': n_sn}
    print(f'{cond_name}: n={n}  (Into_Chromatin={n_ch}, Into_SN={n_sn})  width={w:.2f}"')

max_col_w = max(d['w'] for d in cond_data.values())

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
FRAC_ROW_H   = 0.50    # height per fraction row (Chrom or SN)
BAND_H       = 2 * FRAC_ROW_H      # Chrom + SN
ANNOT_ROW_H  = 0.13    # thin direction strip above each band
SECTION_GAP  = 0.35    # vertical gap between condition bands
N_COND       = 3

LEFT_IN      = 2.60
RIGHT_IN     = 2.80
TOP_PAD      = 0.55
BOT_PAD      = 1.20
CBAR_W_IN    = 2.8
CBAR_H_IN    = 0.14
CBAR_BOT     = 0.46

FIG_W = LEFT_IN + max_col_w + RIGHT_IN
FIG_H = (TOP_PAD
         + N_COND * (ANNOT_ROW_H + BAND_H) + (N_COND - 1) * SECTION_GAP
         + BOT_PAD)

fig      = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
y_cursor = FIG_H - TOP_PAD

COMP_BG = {'Chromatin': '#FBF9F0', 'SN': '#F0F6FB'}

for ci, (cond_name, _, _) in enumerate(CONDITIONS):
    d       = cond_data[cond_name]
    df_cond = d['df']
    n       = d['n']
    col_w   = d['w']          # actual pixel width of this band
    cond_label = COND_LABELS[ci]

    # ── Direction annotation strip ────────────────────────────────────────────
    annot_row = np.where(df_cond['true_dirn'] == 'Into_Chromatin', 1.0, -1.0).reshape(1, -1)
    y_annot_bot = y_cursor - ANNOT_ROW_H
    ax_ann = fig.add_axes([LEFT_IN/FIG_W, y_annot_bot/FIG_H,
                            col_w/FIG_W,   ANNOT_ROW_H/FIG_H])
    ax_ann.imshow(annot_row, aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM,
                  interpolation='nearest')
    ax_ann.set_xticks([]); ax_ann.set_yticks([])
    for sp in ax_ann.spines.values(): sp.set_visible(False)

    y_cursor = y_annot_bot

    # ── Chromatin row ─────────────────────────────────────────────────────────
    ch_mat = df_cond[f'Chrom_FC_{cond_name}'].fillna(0).values.reshape(1, -1)
    ax_ch  = fig.add_axes([LEFT_IN/FIG_W,
                            (y_cursor - FRAC_ROW_H)/FIG_H,
                            col_w/FIG_W,
                            FRAC_ROW_H/FIG_H],
                           facecolor=COMP_BG['Chromatin'])
    im = ax_ch.imshow(ch_mat, aspect='auto', cmap=CMAP, norm=norm,
                      interpolation='nearest')
    ax_ch.set_xticks([]); ax_ch.set_yticks([])
    for sp in ax_ch.spines.values(): sp.set_visible(False)
    # Fraction label (right side, at chromatin row midpoint)
    fig.text((LEFT_IN + col_w + 0.10)/FIG_W,
             (y_cursor - FRAC_ROW_H/2)/FIG_H,
             'Chromatin', ha='left', va='center',
             fontsize=7.5, fontfamily='Arial', color='#111111')

    y_cursor -= FRAC_ROW_H

    # ── SN row ────────────────────────────────────────────────────────────────
    sn_mat = df_cond[f'SN_FC_{cond_name}'].fillna(0).values.reshape(1, -1)
    ax_sn  = fig.add_axes([LEFT_IN/FIG_W,
                            (y_cursor - FRAC_ROW_H)/FIG_H,
                            col_w/FIG_W,
                            FRAC_ROW_H/FIG_H],
                           facecolor=COMP_BG['SN'])
    ax_sn.imshow(sn_mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_sn.set_xticks([]); ax_sn.set_yticks([])
    for sp in ax_sn.spines.values(): sp.set_visible(False)
    fig.text((LEFT_IN + col_w + 0.10)/FIG_W,
             (y_cursor - FRAC_ROW_H/2)/FIG_H,
             'Soluble\nNuclear', ha='left', va='center',
             fontsize=7.5, fontfamily='Arial', color='#111111')

    y_cursor -= FRAC_ROW_H

    # ── Separator between fractions ───────────────────────────────────────────
    sep_y = (y_cursor + FRAC_ROW_H) / FIG_H
    fig.add_artist(plt.Line2D(
        [LEFT_IN/FIG_W, (LEFT_IN + col_w)/FIG_W], [sep_y, sep_y],
        transform=fig.transFigure, color='#CCCCCC', linewidth=0.8, zorder=5))

    # ── Condition label on left ───────────────────────────────────────────────
    mid_y = (y_cursor + FRAC_ROW_H + ANNOT_ROW_H/2 + BAND_H/2) / FIG_H
    fig.text((LEFT_IN - 0.12)/FIG_W, mid_y,
             cond_label, ha='right', va='center',
             fontsize=10, fontweight='bold', fontfamily='Arial', color='#111111')
    # n= breakdown
    fig.text((LEFT_IN - 0.12)/FIG_W, mid_y - 0.058,
             f'n={n}  (↑Chrom={d["n_ch"]}, ↑SN={d["n_sn"]})',
             ha='right', va='top',
             fontsize=7, fontfamily='Arial', color='#555555', style='italic')

    # ── Right-edge tick indicating end of this condition's proteins ───────────
    end_x  = (LEFT_IN + col_w) / FIG_W
    bot_y  = y_cursor / FIG_H
    top_y  = (y_cursor + BAND_H + ANNOT_ROW_H) / FIG_H
    fig.add_artist(plt.Line2D(
        [end_x, end_x], [bot_y, top_y],
        transform=fig.transFigure, color='#AAAAAA', linewidth=0.8,
        linestyle='--', zorder=5))

    y_cursor -= SECTION_GAP

# ── Legend ────────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color='#E8305A', label='Into Chromatin'),
    mpatches.Patch(color='#2B7FD4', label='Into SN'),
]
fig.legend(handles=legend_patches, loc='lower right',
           bbox_to_anchor=(0.98, CBAR_BOT/FIG_H + 0.045),
           fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
           title='True Translocation\n(direction strip)', title_fontsize=7)

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
ax_cb = fig.add_axes([cbar_left, CBAR_BOT/FIG_H, CBAR_W_IN/FIG_W, CBAR_H_IN/FIG_H])
cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
cb.set_ticks([-3,-2,-1,0,1,2,3]); cb.ax.tick_params(labelsize=7)
ax_cb.set_xlabel('Log2 Fold Change vs Naive', fontsize=8, fontfamily='Arial',
                 color='#444444', labelpad=4)

# ── Title & footnote ──────────────────────────────────────────────────────────
fig.text(0.5, (FIG_H - TOP_PAD*0.35)/FIG_H,
         'True Protein Translocation — Per-Condition  [p < 0.10]',
         ha='center', va='top', fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')
fig.text(0.5, 0.008,
         'Each band shows only proteins with true translocation in that condition, '
         'sorted Into_Chromatin → Into_SN  |  '
         'Band width proportional to protein count  |  AW-M-3 excluded  |  Chromatin: Keep+Review only',
         ha='center', va='bottom', fontsize=5.5,
         fontfamily='Arial', color='#888888')

# ── Save ──────────────────────────────────────────────────────────────────────
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
