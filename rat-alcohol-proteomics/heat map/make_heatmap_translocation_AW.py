"""
Translocation heatmap:
- Each ROW = one protein
- Column 1 = Soluble Nuclear AW fold change vs Naïve
- Column 2 = Chromatin AW fold change vs Naïve
- Proteins = union significant in AW in either compartment (n=1279)
- Sorted by Chromatin AW fold change (top = most decreased, bottom = most increased)
- AW-M-3 excluded. Chromatin: Keep+Review filter.
Output: Heatmap_Translocation_AW.pdf
"""

import pandas as pd
import numpy as np
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmap_Translocation_AW.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

CMAP = LinearSegmentedColormap.from_list(
    'bwr',
    ['#2166AC', '#4393C3', '#92C5DE', '#D1E5F0',
     '#FFFFFF',
     '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
    N=512
)

NAIVE_SUFFIXES = ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']
AW_SUFFIXES    = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pvs   = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pvs.append(stats.ttest_ind(n, c, equal_var=False)[1]
                   if len(n) >= 2 and len(c) >= 2 else np.nan)
    pvs   = pd.Series(pvs, index=df.index)
    valid = pvs.notna()
    ranks = pvs[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pvs[valid] * valid.sum() / ranks))
    return fc, corr

# ── Chromatin ──────────────────────────────────────────────────────────────────
raw_c = pd.read_excel(FILE, sheet_name='Chromatin')
df_c  = raw_c[raw_c['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
pref_c    = 'Chrom_'
n_cols_c  = [c for c in df_c.columns if any(pref_c+s in str(c) for s in NAIVE_SUFFIXES)]
aw_cols_c = [c for c in df_c.columns if any(pref_c+s in str(c) for s in AW_SUFFIXES)]
fc_c, corr_c = calc_stats(df_c, n_cols_c, aw_cols_c)
sig_c = (corr_c.notna() & fc_c.notna() & (corr_c > CORR_THRESH) & (fc_c.abs() > FC_THRESH))
# FC map includes ALL proteins (not just significant) so union proteins show a value
chrom_fc_map = dict(zip(df_c['Gene symbol'].astype(str), fc_c.values))
chrom_sig    = set(df_c.loc[sig_c, 'Gene symbol'].astype(str))
print(f'Chromatin AW sig: {len(chrom_sig)}')

# ── Soluble Nuclear ────────────────────────────────────────────────────────────
raw_n = pd.read_excel(FILE, sheet_name='Soluble nuclear')
df_n  = raw_n.copy().reset_index(drop=True)
pref_n    = 'Nuc_'
n_cols_n  = [c for c in df_n.columns if any(pref_n+s in str(c) for s in NAIVE_SUFFIXES)]
aw_cols_n = [c for c in df_n.columns if any(pref_n+s in str(c) for s in AW_SUFFIXES)]
fc_n, corr_n = calc_stats(df_n, n_cols_n, aw_cols_n)
sig_n = (corr_n.notna() & fc_n.notna() & (corr_n > CORR_THRESH) & (fc_n.abs() > FC_THRESH))
# FC map includes ALL proteins so union proteins show a value even if not significant
nuc_fc_map = dict(zip(df_n['Gene symbol'].astype(str), fc_n.values))
nuc_sig    = set(df_n.loc[sig_n, 'Gene symbol'].astype(str))
print(f'Soluble Nuclear AW sig: {len(nuc_sig)}')

# ── Union sig list — only keep proteins detected in BOTH compartments ──────────
# (so both columns are fully filled — no white gaps)
all_chrom_genes = set(df_c['Gene symbol'].astype(str))
all_nuc_genes   = set(df_n['Gene symbol'].astype(str))

union_genes = sorted(
    (chrom_sig | nuc_sig) & all_chrom_genes & all_nuc_genes,
    key=lambda g: chrom_fc_map.get(g, np.nan)
    if not np.isnan(chrom_fc_map.get(g, np.nan)) else np.inf
)
n = len(union_genes)
print(f'Union (detected in both compartments): {n} proteins')

# Matrix: (n_proteins, 2)  col0=SN, col1=Chromatin
sn_col    = np.array([nuc_fc_map.get(g,   np.nan) for g in union_genes])
chrom_col = np.array([chrom_fc_map.get(g, np.nan) for g in union_genes])
mat = np.column_stack([sn_col, chrom_col])   # (n_proteins, 2)

# ── Figure — tall, narrow (2 columns, many rows) ───────────────────────────────
ROW_H_IN  = 0.013   # inches per protein row (same as individual heatmap scripts)
COL_W_IN  = 0.8     # width per column
LEFT_PAD  = 1.80
RIGHT_PAD = 0.40
TOP_PAD   = 0.70
BOT_PAD   = 0.80

hm_h = n * ROW_H_IN
hm_w = 2 * COL_W_IN

FIG_W = LEFT_PAD + hm_w + RIGHT_PAD
FIG_H = TOP_PAD + hm_h + BOT_PAD

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')

vmax = 2.0
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

def ffy(y): return y / FIG_H
def ffx(x): return x / FIG_W

band_bot = BOT_PAD
ax = fig.add_axes([ffx(LEFT_PAD), ffy(band_bot), ffx(hm_w), ffy(hm_h)])
im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')

# White divider between the 2 columns
ax.axvline(0.5, color='white', linewidth=2.0, zorder=3)

ax.set_yticks([])
ax.set_xticks([0, 1])
ax.set_xticklabels(['Soluble\nNuclear', 'Chromatin'],
                   fontsize=11, fontweight='bold', fontfamily='Arial')
ax.xaxis.set_tick_params(length=0)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()

for sp in ax.spines.values(): sp.set_visible(False)

# Title
fig.text(ffx(LEFT_PAD + hm_w / 2),
         ffy(band_bot + hm_h + 0.40),
         f'Acute Withdrawal vs Naïve  (n = {n} proteins)',
         ha='center', va='bottom',
         fontsize=12, fontweight='bold', fontfamily='Arial', color='#111111')

# Colorbar
cbar_bot = 0.15
ax_cbar = fig.add_axes([ffx(LEFT_PAD - 0.3), ffy(cbar_bot), ffx(hm_w + 0.6), ffy(0.14)])
cb = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
cb.set_label('')
cb.ax.tick_params(labelsize=8)
cb.set_ticks([-2, -1, 0, 1, 2])
ax_cbar.text(0.5, -0.9, 'Log\u2082 Fold Change vs Na\u00efve',
             transform=ax_cbar.transAxes,
             ha='center', va='top', fontsize=8.5,
             fontfamily='Arial', color='#444444')

pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
pdf.close()
plt.close(fig)
print(f'Saved: {OUT}')
