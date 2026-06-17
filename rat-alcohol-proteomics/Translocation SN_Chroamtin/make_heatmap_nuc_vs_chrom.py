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

FILE     = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_PDF  = 'Heatmap_Nuc_vs_Chrom_AW.pdf'
OUT_XLSX = 'Nuc_vs_Chrom_AW_stats.xlsx'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
ROW_H_IN    = 0.013

# Red-white-blue colormap
CMAP = LinearSegmentedColormap.from_list(
    'rwb',
    ['#2166AC', '#74ADD1', '#ABD9E9', '#FFFFFF', '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
    N=256
)

# Replicate suffixes per condition
COND_SUFFIXES = {
    'Naive':                 ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
    'Protracted Abstinence': ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
AW_SUFFIXES = COND_SUFFIXES['Acute Withdrawal']

COMP_COLORS = {
    'Soluble nuclear': '#E8DAEF',
    'Chromatin':       '#D5F5E3',
}

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
    p_vals    = pd.Series(p_vals, index=df.index)
    valid     = p_vals.notna()
    ranks     = p_vals[valid].rank(ascending=False)
    corrected = pd.Series(np.nan, index=df.index)
    corrected[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corrected

def load_compartment(excel_name, prefix):
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    df  = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    naive_cols = [c for c in df.columns if prefix + 'N-' in str(c)]
    aw_cols    = [c for c in df.columns if any(prefix + s in str(c) for s in AW_SUFFIXES)]
    # All condition cols
    cond_cols_all = {}
    for cname, suffixes in COND_SUFFIXES.items():
        cond_cols_all[cname] = [c for c in df.columns
                                if any(prefix + s in str(c) for s in suffixes)]
    return df, naive_cols, aw_cols, cond_cols_all

def get_values(lookup, gene, cols):
    if gene not in lookup.index:
        return np.array([np.nan] * len(cols))
    row = lookup.loc[gene]
    return np.array([float(pd.to_numeric(row[c], errors='coerce'))
                     if c in row.index else np.nan for c in cols])

# ── Load both compartments ──
nuc_df,   nuc_naive,   nuc_aw,   nuc_cond_cols   = load_compartment('Soluble nuclear', 'Nuc_')
chrom_df, chrom_naive, chrom_aw, chrom_cond_cols = load_compartment('Chromatin',       'Chrom_')

# Significant in AW vs Naive
nuc_fc,   nuc_corr   = calc_stats(nuc_df,   nuc_naive,   nuc_aw)
chrom_fc, chrom_corr = calc_stats(chrom_df, chrom_naive, chrom_aw)

nuc_sig_mask   = (nuc_corr   > CORR_THRESH) & (nuc_fc.abs()   > FC_THRESH)
chrom_sig_mask = (chrom_corr > CORR_THRESH) & (chrom_fc.abs() > FC_THRESH)

nuc_sig_genes   = set(nuc_df.loc[nuc_sig_mask,    'Gene symbol'].astype(str))
chrom_sig_genes = set(chrom_df.loc[chrom_sig_mask, 'Gene symbol'].astype(str))
union_genes     = nuc_sig_genes | chrom_sig_genes

print(f'Nuc AW sig:   {len(nuc_sig_genes)}')
print(f'Chrom AW sig: {len(chrom_sig_genes)}')
print(f'Union:        {len(union_genes)}')

# Index compartments by gene symbol
nuc_lookup   = (nuc_df.drop_duplicates('Gene symbol')
                .set_index(nuc_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))
chrom_lookup = (chrom_df.drop_duplicates('Gene symbol')
                .set_index(chrom_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))

# ── Build heatmap matrix (AW only, 10 cols) ──
all_col_names = nuc_aw + chrom_aw
gene_list = sorted(union_genes)
n_genes   = len(gene_list)
n_cols    = len(all_col_names)

mat = np.full((n_genes, n_cols), np.nan)
for row_i, gene in enumerate(gene_list):
    nv = get_values(nuc_lookup,   gene, nuc_aw)
    cv = get_values(chrom_lookup, gene, chrom_aw)
    mat[row_i, :len(nuc_aw)]  = nv
    mat[row_i, len(nuc_aw):]  = cv

# Z-score each row
row_mean = np.nanmean(mat, axis=1, keepdims=True)
row_std  = np.nanstd(mat,  axis=1, keepdims=True)
row_std[row_std == 0] = 1
mat_z = (mat - row_mean) / row_std

mat_fill = np.where(np.isnan(mat_z), 0, mat_z)
if n_genes > 1:
    Z = linkage(mat_fill, method='ward')
    order = leaves_list(Z)
else:
    order = [0]; Z = None

genes_ordered = [gene_list[i] for i in order]
mat_ordered   = mat_z[order, :]
print(f'Matrix shape: {mat_ordered.shape}')

# ── Statistics: per-protein Chrom−Nuc FC per condition ──
COND_ORDER = ['Naive', 'Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

stat_rows = []
fc_by_cond = {c: [] for c in COND_ORDER}   # collect all per-protein FCs per condition

for gene in gene_list:
    row = {'Gene': gene,
           'Sig_in': ('Both'     if gene in nuc_sig_genes and gene in chrom_sig_genes else
                      'Nuc only' if gene in nuc_sig_genes else 'Chrom only')}
    for cname in COND_ORDER:
        nuc_vals   = get_values(nuc_lookup,   gene, nuc_cond_cols[cname])
        chrom_vals = get_values(chrom_lookup, gene, chrom_cond_cols[cname])
        nuc_mean   = np.nanmean(nuc_vals)
        chrom_mean = np.nanmean(chrom_vals)
        fc = chrom_mean - nuc_mean   # log-space Chrom minus Nuc
        row[f'FC_Chrom_vs_Nuc_{cname}'] = round(fc, 4) if np.isfinite(fc) else np.nan
        if np.isfinite(fc):
            fc_by_cond[cname].append(fc)
    stat_rows.append(row)

stat_df = pd.DataFrame(stat_rows)

# One-way ANOVA: do mean Chrom/Nuc FCs differ across the 4 conditions?
# Each group = all per-protein FCs for that condition
anova_groups = [np.array(fc_by_cond[c]) for c in COND_ORDER]
f_stat, anova_pval = stats.f_oneway(*anova_groups)

# Summary stats per condition
summary_rows = []
for cname in COND_ORDER:
    vals = np.array(fc_by_cond[cname])
    summary_rows.append({
        'Condition':  cname,
        'N_proteins': len(vals),
        'Mean_FC_Chrom_vs_Nuc': round(np.mean(vals), 4),
        'SEM':        round(stats.sem(vals), 4),
        'SD':         round(np.std(vals, ddof=1), 4),
    })
summary_df = pd.DataFrame(summary_rows)
summary_df['ANOVA_F'] = round(f_stat, 4)
summary_df['ANOVA_pval'] = anova_pval
print(f'\nOne-way ANOVA: F={f_stat:.4f}, p={anova_pval:.4e}')
print(summary_df[['Condition','Mean_FC_Chrom_vs_Nuc','SEM','ANOVA_pval']].to_string(index=False))

# Save Excel with two sheets
with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    summary_df.to_excel(writer, sheet_name='Summary', index=False)
    stat_df.to_excel(writer, sheet_name='Per_protein_FC', index=False)
print(f'Saved stats: {OUT_XLSX}')

# ── Bar plot: mean Chrom/Nuc FC per condition ──
fig_bar, ax_bar = plt.subplots(figsize=(7, 5))
fig_bar.patch.set_facecolor('white')
cond_labels  = ['Naïve', 'Intoxication', 'Acute\nWithdrawal', 'Protracted\nAbstinence']
bar_colors   = ['#D6EAF8', '#FDEBD0', '#D5F5E3', '#E8DAEF']
means = [summary_df.loc[i, 'Mean_FC_Chrom_vs_Nuc'] for i in range(4)]
sems  = [summary_df.loc[i, 'SEM']                  for i in range(4)]

bars = ax_bar.bar(range(4), means, yerr=sems, color=bar_colors,
                  edgecolor='#333333', linewidth=1.5,
                  error_kw=dict(ecolor='#333333', lw=2, capsize=6, capthick=2))
ax_bar.axhline(0, color='#888888', linewidth=1.0, linestyle='--')
ax_bar.set_xticks(range(4))
ax_bar.set_xticklabels(cond_labels, fontsize=12)
ax_bar.set_ylabel('Mean FC  (Chromatin − Nuclear, log₂)', fontsize=12)
ax_bar.set_title('Average Chromatin/Nuclear protein ratio\nby condition', fontsize=13, fontweight='bold')
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)

# ANOVA p-value annotation
p_txt = f'One-way ANOVA  p = {anova_pval:.2e}'
ax_bar.text(0.98, 0.97, p_txt, transform=ax_bar.transAxes,
            ha='right', va='top', fontsize=11, color='#333333',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='#CCCCCC', alpha=0.9))

plt.tight_layout()
fig_bar.savefig('Chrom_vs_Nuc_mean_FC_barplot.pdf', dpi=150, bbox_inches='tight')
plt.close(fig_bar)
print('Saved: Chrom_vs_Nuc_mean_FC_barplot.pdf')

# ── Heatmap ──
COL_W    = 0.26
LABEL_W  = 1.5
DENDRO_W = 0.5
HEADER_H = 0.9
TITLE_H  = 0.6
BOTTOM   = 0.4

hm_h  = n_genes * ROW_H_IN
fig_h = hm_h + HEADER_H + TITLE_H + BOTTOM + 0.6
fig_w = DENDRO_W + LABEL_W + n_cols * COL_W + 0.8

fig = plt.figure(figsize=(fig_w, fig_h))
fig.patch.set_facecolor('white')

vabs = float(np.nanpercentile(np.abs(mat_ordered[np.isfinite(mat_ordered)]), 96))
vabs = max(vabs, 0.5)

x0          = (DENDRO_W + LABEL_W + 0.3) / fig_w
heat_w_f    = (n_cols * COL_W) / fig_w
heat_top    = (fig_h - TITLE_H - HEADER_H - 0.05) / fig_h
heat_bot    = heat_top - hm_h / fig_h
label_left  = (DENDRO_W + 0.3) / fig_w
label_w_f   = LABEL_W / fig_w
dendro_left = 0.3 / fig_w
dendro_w_f  = DENDRO_W / fig_w

ax_heat   = fig.add_axes([x0,          heat_bot, heat_w_f,          hm_h / fig_h])
ax_labels = fig.add_axes([label_left,  heat_bot, label_w_f,         hm_h / fig_h])
ax_dendro = fig.add_axes([dendro_left, heat_bot, dendro_w_f * 0.85, hm_h / fig_h])
ax_header = fig.add_axes([x0,          heat_top, heat_w_f, HEADER_H * 0.65 / fig_h])

norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
im   = ax_heat.imshow(mat_ordered, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
ax_heat.set_yticks([])
ax_heat.set_xticks([])
for sp in ax_heat.spines.values():
    sp.set_visible(False)

ax_heat.axvline(len(nuc_aw) - 0.5, color='white', linewidth=3.0)

for j, col in enumerate(all_col_names):
    short = str(col).split('_')[-1] if '_' in str(col) else str(col)
    ax_heat.text(j, n_genes + 0.5, short,
                 ha='center', va='top', fontsize=5.5, fontfamily='Arial',
                 transform=ax_heat.transData, rotation=55, color='#444444')

if Z is not None and n_genes > 1:
    dendrogram(Z, ax=ax_dendro, orientation='left',
               color_threshold=0, above_threshold_color='#777777',
               no_labels=True, count_sort='ascending')
ax_dendro.set_xlim(ax_dendro.get_xlim()[::-1])
ax_dendro.axis('off')

font_size = max(2.5, min(7, 350 / n_genes))
ax_labels.set_xlim(0, 1)
ax_labels.set_ylim(-0.5, n_genes - 0.5)
ax_labels.invert_yaxis()
for i, gene in enumerate(genes_ordered):
    ax_labels.text(0.95, i, gene, ha='right', va='center',
                   fontsize=font_size, fontfamily='Arial', color='#222222')
ax_labels.axis('off')

blocks = [
    ('Soluble Nuclear\n(AW replicates)', COMP_COLORS['Soluble nuclear'], 0,           len(nuc_aw)),
    ('Chromatin\n(AW replicates)',        COMP_COLORS['Chromatin'],       len(nuc_aw), len(nuc_aw) + len(chrom_aw)),
]
ax_header.set_xlim(0, n_cols)
ax_header.set_ylim(0, 1)
for name, color, start, end in blocks:
    nc = end - start
    ax_header.add_patch(mpatches.FancyBboxPatch(
        (start + 0.08, 0.1), nc - 0.16, 0.88,
        boxstyle='round,pad=0.02', facecolor=color,
        edgecolor='white', linewidth=2.0))
    ax_header.text(start + nc / 2, 0.54, name,
                   ha='center', va='center', fontsize=9,
                   fontweight='bold', fontfamily='Arial')
ax_header.axis('off')

cbar_left = x0 + heat_w_f + 0.01
ax_cbar   = fig.add_axes([cbar_left, heat_bot + hm_h * 0.3 / fig_h,
                           0.015, hm_h * 0.4 / fig_h])
cb = fig.colorbar(im, cax=ax_cbar)
cb.set_label('Z-score', fontsize=8, fontfamily='Arial')
cb.ax.tick_params(labelsize=7)

fig.text(0.5, (fig_h - TITLE_H * 0.35) / fig_h,
         f'Soluble Nuclear vs Chromatin — Acute Withdrawal replicates\n'
         f'Union of AW-significant proteins (corrected p > {CORR_THRESH}, |FC| > {FC_THRESH})   n={n_genes}',
         ha='center', va='center', fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#222222')

plt.savefig(OUT_PDF, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {OUT_PDF}')
