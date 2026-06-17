import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE    = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_PDF = 'Heatmap_Nuc_Chrom_AllConditions.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
ROW_H_IN    = 0.013

# Red-white-blue colormap
CMAP = LinearSegmentedColormap.from_list(
    'rwb',
    ['#2166AC', '#74ADD1', '#ABD9E9', '#FFFFFF', '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
    N=256
)

COND_SUFFIXES = {
    'Naive':                 ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
    'Protracted Abstinence': ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
COND_ORDER  = ['Naive', 'Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
COND_COLORS = ['#D6EAF8', '#FDEBD0', '#D5F5E3', '#E8DAEF']
AW_SUFFIXES = COND_SUFFIXES['Acute Withdrawal']

# Compartment header colors
COMP_COLORS = {'Chromatin': '#A9DFBF', 'Soluble nuclear': '#D7BDE2'}

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        if len(n) >= 2 and len(c) >= 2:
            _, p = stats.ttest_ind(n, c, equal_var=False)
            p_vals.append(p)
        else:
            p_vals.append(np.nan)
    p_vals  = pd.Series(p_vals, index=df.index)
    valid   = p_vals.notna()
    ranks   = p_vals[valid].rank(ascending=False)
    corr    = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corr

def load_compartment(excel_name, prefix):
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    df  = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    naive_cols = [c for c in df.columns if prefix + 'N-' in str(c)]
    aw_cols    = [c for c in df.columns if any(prefix + s in str(c) for s in AW_SUFFIXES)]
    # All rep cols in condition order
    all_rep_cols = []
    cond_col_groups = []
    for cname in COND_ORDER:
        gcols = [c for c in df.columns
                 if any(prefix + s in str(c) for s in COND_SUFFIXES[cname])]
        cond_col_groups.append((cname, gcols))
        all_rep_cols.extend(gcols)
    return df, naive_cols, aw_cols, all_rep_cols, cond_col_groups

def get_values(lookup, gene, cols):
    if gene not in lookup.index:
        return np.full(len(cols), np.nan)
    row = lookup.loc[gene]
    return np.array([float(pd.to_numeric(row[c], errors='coerce'))
                     if c in row.index else np.nan for c in cols])

# ── Load both compartments ──
nuc_df,   nuc_naive,   nuc_aw,   nuc_rep_cols,   nuc_groups   = load_compartment('Soluble nuclear', 'Nuc_')
chrom_df, chrom_naive, chrom_aw, chrom_rep_cols, chrom_groups = load_compartment('Chromatin',       'Chrom_')

# Union of AW-significant proteins
nuc_fc,   nuc_corr   = calc_stats(nuc_df,   nuc_naive, nuc_aw)
chrom_fc, chrom_corr = calc_stats(chrom_df, chrom_naive, chrom_aw)
nuc_sig   = set(nuc_df.loc[(nuc_corr   > CORR_THRESH) & (nuc_fc.abs()   > FC_THRESH), 'Gene symbol'].astype(str))
chrom_sig = set(chrom_df.loc[(chrom_corr > CORR_THRESH) & (chrom_fc.abs() > FC_THRESH), 'Gene symbol'].astype(str))
union_genes = sorted(nuc_sig | chrom_sig)
n_genes = len(union_genes)
print(f'Union: {n_genes} proteins')

nuc_lookup   = (nuc_df.drop_duplicates('Gene symbol')
                .set_index(nuc_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))
chrom_lookup = (chrom_df.drop_duplicates('Gene symbol')
                .set_index(chrom_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))

# ── Build matrix: cols = [Chrom x 20, Nuc x 20] ──
# Column layout: Chrom (N|I|AW|PA) then Nuc (N|I|AW|PA)
all_cols     = chrom_rep_cols + nuc_rep_cols
n_chrom_cols = len(chrom_rep_cols)
n_nuc_cols   = len(nuc_rep_cols)
n_cols       = len(all_cols)
print(f'Columns: {n_chrom_cols} Chromatin + {n_nuc_cols} Nuclear = {n_cols} total')

mat = np.full((n_genes, n_cols), np.nan)
for row_i, gene in enumerate(union_genes):
    chrom_vals = get_values(chrom_lookup, gene, chrom_rep_cols)
    nuc_vals   = get_values(nuc_lookup,   gene, nuc_rep_cols)
    mat[row_i, :n_chrom_cols] = chrom_vals
    mat[row_i, n_chrom_cols:] = nuc_vals

# Z-score each row across all 40 columns
row_mean = np.nanmean(mat, axis=1, keepdims=True)
row_std  = np.nanstd(mat,  axis=1, keepdims=True)
row_std[row_std == 0] = 1
mat_z = (mat - row_mean) / row_std

# Cluster rows
mat_fill = np.where(np.isnan(mat_z), 0, mat_z)
if n_genes > 1:
    Z = linkage(mat_fill, method='ward')
    order = leaves_list(Z)
else:
    order = [0]; Z = None

genes_ordered = [union_genes[i] for i in order]
mat_ordered   = mat_z[order, :]

# ── Figure layout ──
COL_W    = 0.15    # narrower since more columns
LABEL_W  = 1.5
DENDRO_W = 0.5
# Two-tier header: top = compartment, bottom = condition
HEADER1_H = 0.45   # compartment bar
HEADER2_H = 0.55   # condition bars
TITLE_H   = 0.6
BOTTOM    = 0.4

hm_h  = n_genes * ROW_H_IN
fig_h = hm_h + HEADER1_H + HEADER2_H + TITLE_H + BOTTOM + 0.5
fig_w = DENDRO_W + LABEL_W + n_cols * COL_W + 0.8

fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor('white')

vabs = float(np.nanpercentile(np.abs(mat_ordered[np.isfinite(mat_ordered)]), 96))
vabs = max(vabs, 0.5)

x0          = (DENDRO_W + LABEL_W + 0.3) / fig_w
heat_w_f    = (n_cols * COL_W) / fig_w
heat_top    = (fig_h - TITLE_H - HEADER1_H - HEADER2_H - 0.05) / fig_h
heat_bot    = heat_top - hm_h / fig_h
label_left  = (DENDRO_W + 0.3) / fig_w
label_w_f   = LABEL_W / fig_w
dendro_left = 0.3 / fig_w

ax_heat    = fig.add_axes([x0,           heat_bot, heat_w_f,          hm_h / fig_h])
ax_labels  = fig.add_axes([label_left,   heat_bot, label_w_f,         hm_h / fig_h])
ax_dendro  = fig.add_axes([dendro_left,  heat_bot, DENDRO_W*0.85/fig_w, hm_h / fig_h])
# Condition header (lower tier)
ax_cond    = fig.add_axes([x0, heat_top,                    heat_w_f, HEADER2_H * 0.7 / fig_h])
# Compartment header (upper tier)
ax_comp    = fig.add_axes([x0, heat_top + HEADER2_H/fig_h, heat_w_f, HEADER1_H * 0.8 / fig_h])

# ── Heatmap ──
norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
im   = ax_heat.imshow(mat_ordered, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
ax_heat.set_yticks([]); ax_heat.set_xticks([])
for sp in ax_heat.spines.values(): sp.set_visible(False)

# White divider between Chrom and Nuc blocks
ax_heat.axvline(n_chrom_cols - 0.5, color='white', linewidth=4.0)

# Lighter dividers between conditions within each compartment
col_cursor = 0
for _, gcols in chrom_groups:
    col_cursor += len(gcols)
    if col_cursor < n_chrom_cols:
        ax_heat.axvline(col_cursor - 0.5, color='white', linewidth=1.5)
col_cursor = n_chrom_cols
for _, gcols in nuc_groups:
    col_cursor += len(gcols)
    if col_cursor < n_cols:
        ax_heat.axvline(col_cursor - 0.5, color='white', linewidth=1.5)

# Column rep labels
for j, col in enumerate(all_cols):
    short = str(col).split('_')[-1] if '_' in str(col) else str(col)
    ax_heat.text(j, n_genes + 0.5, short,
                 ha='center', va='top', fontsize=4.5, fontfamily='Arial',
                 transform=ax_heat.transData, rotation=60, color='#444444')

# ── Dendrogram ──
if Z is not None and n_genes > 1:
    dendrogram(Z, ax=ax_dendro, orientation='left',
               color_threshold=0, above_threshold_color='#777777',
               no_labels=True, count_sort='ascending')
ax_dendro.set_xlim(ax_dendro.get_xlim()[::-1])
ax_dendro.axis('off')

# ── Gene labels ──
font_size = max(2.5, min(7, 350 / n_genes))
ax_labels.set_xlim(0, 1); ax_labels.set_ylim(-0.5, n_genes - 0.5)
ax_labels.invert_yaxis()
for i, gene in enumerate(genes_ordered):
    ax_labels.text(0.95, i, gene, ha='right', va='center',
                   fontsize=font_size, fontfamily='Arial', color='#222222')
ax_labels.axis('off')

# ── Condition header (lower tier) ──
ax_cond.set_xlim(0, n_cols); ax_cond.set_ylim(0, 1)
# Chrom conditions
col_cursor = 0
for (cname, gcols), ccolor in zip(chrom_groups, COND_COLORS):
    nc = len(gcols)
    ax_cond.add_patch(mpatches.FancyBboxPatch(
        (col_cursor + 0.06, 0.08), nc - 0.12, 0.84,
        boxstyle='round,pad=0.02', facecolor=ccolor, edgecolor='white', linewidth=1.2))
    short = cname.replace('Intoxication','Intox').replace('Acute Withdrawal','AW').replace('Protracted Abstinence','PA')
    ax_cond.text(col_cursor + nc/2, 0.52, short,
                 ha='center', va='center', fontsize=6.5, fontweight='bold', fontfamily='Arial')
    col_cursor += nc
# Nuc conditions
for (cname, gcols), ccolor in zip(nuc_groups, COND_COLORS):
    nc = len(gcols)
    ax_cond.add_patch(mpatches.FancyBboxPatch(
        (col_cursor + 0.06, 0.08), nc - 0.12, 0.84,
        boxstyle='round,pad=0.02', facecolor=ccolor, edgecolor='white', linewidth=1.2))
    short = cname.replace('Intoxication','Intox').replace('Acute Withdrawal','AW').replace('Protracted Abstinence','PA')
    ax_cond.text(col_cursor + nc/2, 0.52, short,
                 ha='center', va='center', fontsize=6.5, fontweight='bold', fontfamily='Arial')
    col_cursor += nc
ax_cond.axis('off')

# ── Compartment header (upper tier) ──
ax_comp.set_xlim(0, n_cols); ax_comp.set_ylim(0, 1)
for start, width, comp_name in [
    (0,            n_chrom_cols, 'Chromatin'),
    (n_chrom_cols, n_nuc_cols,   'Soluble Nuclear'),
]:
    ax_comp.add_patch(mpatches.FancyBboxPatch(
        (start + 0.1, 0.06), width - 0.2, 0.88,
        boxstyle='round,pad=0.02',
        facecolor=COMP_COLORS[comp_name if comp_name == 'Chromatin' else 'Soluble nuclear'],
        edgecolor='white', linewidth=2.0))
    ax_comp.text(start + width/2, 0.52, comp_name,
                 ha='center', va='center', fontsize=9,
                 fontweight='bold', fontfamily='Arial')
ax_comp.axis('off')

# ── Colorbar ──
cbar_left = x0 + heat_w_f + 0.01
ax_cbar   = fig.add_axes([cbar_left, heat_bot + hm_h*0.3/fig_h, 0.014, hm_h*0.4/fig_h])
cb = fig.colorbar(im, cax=ax_cbar)
cb.set_label('Z-score', fontsize=8, fontfamily='Arial')
cb.ax.tick_params(labelsize=7)

# ── Title ──
fig.text(0.5, (fig_h - TITLE_H*0.35) / fig_h,
         f'Chromatin vs Soluble Nuclear — all conditions  |  '
         f'AW-significant proteins (corrected p > {CORR_THRESH}, |FC| > {FC_THRESH})   n={n_genes}',
         ha='center', va='center', fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#222222')

plt.savefig(OUT_PDF, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {OUT_PDF}')
