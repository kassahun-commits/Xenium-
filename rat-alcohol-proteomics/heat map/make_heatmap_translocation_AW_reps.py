"""
Translocation heatmap — individual replicates version.
Same protein list and layout as Heatmap_Translocation_AW.pdf but:
- 8 columns: SN AW reps (x4) | Chromatin AW reps (x4)
- Color = Z-score of LFQ intensity per protein (z-scored across all 8 AW reps)
- Thick white divider between SN and Chromatin groups
- AW-M-3 excluded. Chromatin: Keep+Review filter.
- Sorted by Chromatin AW fold change (same order as fold change version)
Output: Heatmap_Translocation_AW_Reps.pdf
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
OUT  = 'Heatmap_Translocation_AW_Reps.pdf'

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
AW_SUFFIXES    = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']   # AW-M-3 excluded

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
chrom_fc_map   = dict(zip(df_c['Gene symbol'].astype(str), fc_c.values))
chrom_sig      = set(df_c.loc[sig_c, 'Gene symbol'].astype(str))
chrom_rep_data = df_c.set_index(df_c['Gene symbol'].astype(str))[aw_cols_c].apply(pd.to_numeric, errors='coerce')
print(f'Chromatin AW sig: {len(chrom_sig)}')

# ── Soluble Nuclear ────────────────────────────────────────────────────────────
raw_n = pd.read_excel(FILE, sheet_name='Soluble nuclear')
df_n  = raw_n.copy().reset_index(drop=True)
pref_n    = 'Nuc_'
n_cols_n  = [c for c in df_n.columns if any(pref_n+s in str(c) for s in NAIVE_SUFFIXES)]
aw_cols_n = [c for c in df_n.columns if any(pref_n+s in str(c) for s in AW_SUFFIXES)]
fc_n, corr_n = calc_stats(df_n, n_cols_n, aw_cols_n)
sig_n = (corr_n.notna() & fc_n.notna() & (corr_n > CORR_THRESH) & (fc_n.abs() > FC_THRESH))
nuc_fc_map   = dict(zip(df_n['Gene symbol'].astype(str), fc_n.values))
nuc_sig      = set(df_n.loc[sig_n, 'Gene symbol'].astype(str))
nuc_rep_data = df_n.set_index(df_n['Gene symbol'].astype(str))[aw_cols_n].apply(pd.to_numeric, errors='coerce')
print(f'Soluble Nuclear AW sig: {len(nuc_sig)}')

# ── Union — only proteins detected in both compartments ───────────────────────
all_chrom = set(df_c['Gene symbol'].astype(str))
all_nuc   = set(df_n['Gene symbol'].astype(str))
union_genes = sorted(
    (chrom_sig | nuc_sig) & all_chrom & all_nuc,
    key=lambda g: chrom_fc_map.get(g, np.nan)
    if not np.isnan(chrom_fc_map.get(g, np.nan)) else np.inf
)
n = len(union_genes)
print(f'Union (detected in both): {n} proteins')

# ── Build matrix (n_proteins, 8) ──────────────────────────────────────────────
def get_rep_rows(gene_list, rep_data):
    rows = []
    for g in gene_list:
        if g in rep_data.index:
            vals = rep_data.loc[g].values.astype(float)
            if vals.ndim > 1:
                vals = vals[0]   # handle duplicate gene symbols
        else:
            vals = np.full(len(AW_SUFFIXES), np.nan)
        rows.append(vals)
    return np.array(rows)   # (n_proteins, 4)

mat_nuc   = get_rep_rows(union_genes, nuc_rep_data)    # (n, 4)
mat_chrom = get_rep_rows(union_genes, chrom_rep_data)  # (n, 4)

# Concatenate: SN reps first, then Chromatin reps → (n, 8)
mat_all = np.hstack([mat_nuc, mat_chrom])

# Z-score each protein across all 8 columns
row_mean = np.nanmean(mat_all, axis=1, keepdims=True)
row_std  = np.nanstd(mat_all,  axis=1, keepdims=True)
row_std[row_std == 0] = 1
mat_z = (mat_all - row_mean) / row_std   # (n, 8)

# ── Figure ─────────────────────────────────────────────────────────────────────
ROW_H_IN  = 0.013
COL_W_IN  = 0.40   # wider than fold change version since only 8 cols
LEFT_PAD  = 1.80
RIGHT_PAD = 2.80   # room for column labels on top
TOP_PAD   = 0.90
BOT_PAD   = 0.80

hm_h = n * ROW_H_IN
hm_w = 8 * COL_W_IN

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
im = ax.imshow(mat_z, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')

# Thick white divider between SN (cols 0-3) and Chromatin (cols 4-7)
ax.axvline(3.5, color='white', linewidth=3.0, zorder=3)
# Thin dividers between individual reps
for x in [0.5, 1.5, 2.5, 4.5, 5.5, 6.5]:
    ax.axvline(x, color='white', linewidth=0.6, zorder=2)

# X-axis tick labels = replicate names
ax.set_xticks(range(8))
ax.set_xticklabels(
    [s for s in AW_SUFFIXES] + [s for s in AW_SUFFIXES],
    fontsize=7, fontfamily='Arial', rotation=45, ha='left'
)
ax.xaxis.set_tick_params(length=0)
ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()
ax.set_yticks([])
for sp in ax.spines.values(): sp.set_visible(False)

# Compartment group labels above columns
sn_cx    = ffx(LEFT_PAD + 2 * COL_W_IN)
chrom_cx = ffx(LEFT_PAD + 6 * COL_W_IN)
label_y  = ffy(band_bot + hm_h + 0.55)
fig.text(sn_cx,    label_y, 'Soluble Nuclear',
         ha='center', va='bottom', fontsize=11, fontweight='bold',
         fontfamily='Arial', color='#111111')
fig.text(chrom_cx, label_y, 'Chromatin',
         ha='center', va='bottom', fontsize=11, fontweight='bold',
         fontfamily='Arial', color='#111111')

# Title
fig.text(ffx(LEFT_PAD + hm_w / 2),
         ffy(band_bot + hm_h + 0.78),
         f'Acute Withdrawal — Individual Replicates  (n = {n} proteins)',
         ha='center', va='bottom',
         fontsize=12, fontweight='bold', fontfamily='Arial', color='#111111')

# Colorbar
ax_cbar = fig.add_axes([ffx(LEFT_PAD + hm_w / 2 - 1.25), ffy(0.18), ffx(2.5), ffy(0.14)])
cb = fig.colorbar(im, cax=ax_cbar, orientation='horizontal')
cb.set_label('')
cb.ax.tick_params(labelsize=8)
cb.set_ticks([-2, -1, 0, 1, 2])
ax_cbar.text(0.5, -0.9, 'Z-score (LFQ intensity)',
             transform=ax_cbar.transAxes,
             ha='center', va='top', fontsize=8.5,
             fontfamily='Arial', color='#444444')

fig.text(0.5, ffy(0.04),
         'AW-M-3 excluded. Sorted by Chromatin AW fold change.',
         ha='center', va='bottom', fontsize=7.5,
         color='#888888', style='italic', fontfamily='Arial')

pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
pdf.close()
plt.close(fig)
print(f'Saved: {OUT}')
