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

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmaps_v3.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
ROW_H_IN    = 0.013   # inches per protein row

# Yellow-white-blue colormap (matches paper Fig 4 style)
CMAP = LinearSegmentedColormap.from_list(
    'ywb',
    ['#2166AC', '#74ADD1', '#ABD9E9', '#FFFFFF', '#FEE090', '#FDAE61', '#F46D43'],
    N=256
)

# Condition group header colors (pastel)
GROUP_META = [
    ('Naïve',                  '#D6EAF8', ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']),
    ('Intoxication',           '#FDEBD0', ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2']),
    ('Acute Withdrawal',       '#D5F5E3', ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3']),
    ('Protracted Abstinence',  '#E8DAEF', ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2']),
]

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_',   'keep_review'),
]

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

# ── Load data and compute union significant sets ──
compartment_data = []
for label, excel_name, prefix, fmode in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    if fmode == 'keep_review':
        df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    # Build column lists per group
    groups = []
    all_rep_cols = []
    naive_cols = []
    for gname, gcol, suffixes in GROUP_META:
        gcols = [c for c in df.columns if any(prefix + s in str(c) for s in suffixes)]
        groups.append((gname, gcol, gcols))
        all_rep_cols.extend(gcols)
        if gname == 'Naïve':
            naive_cols = gcols

    # Union of significant proteins across the 3 non-naive conditions
    sig_mask = pd.Series(False, index=df.index)
    for gname, gcol, cond_cols in groups:
        if gname == 'Naïve' or not cond_cols:
            continue
        fc, corr = calc_stats(df, naive_cols, cond_cols)
        sig = corr.notna() & fc.notna() & (corr > CORR_THRESH) & (fc.abs() > FC_THRESH)
        sig_mask = sig_mask | sig

    df_sig = df[sig_mask].reset_index(drop=True)
    n = len(df_sig)
    print(f'{label}: {n} proteins in union')

    if n == 0:
        continue

    mat = df_sig[all_rep_cols].apply(pd.to_numeric, errors='coerce').values.astype(float)
    row_mean = np.nanmean(mat, axis=1, keepdims=True)
    row_std  = np.nanstd(mat,  axis=1, keepdims=True)
    row_std[row_std == 0] = 1
    mat_z = (mat - row_mean) / row_std

    mat_fill = np.where(np.isnan(mat_z), 0, mat_z)
    if n > 1:
        Z = linkage(mat_fill, method='ward')
        order = leaves_list(Z)
    else:
        order = [0]
        Z = None

    genes = df_sig['Gene symbol'].astype(str).iloc[order].tolist()
    mat_ordered = mat_z[order, :]

    compartment_data.append({
        'label':       label,
        'n':           n,
        'mat':         mat_ordered,
        'genes':       genes,
        'groups':      groups,
        'all_rep_cols':all_rep_cols,
        'linkage':     Z,
    })

# ── Figure layout: all 4 side by side, heights proportional ──
N_COMP   = len(compartment_data)
COL_W    = 0.22          # inches per replicate column
LABEL_W  = 1.5           # gene label column width
DENDRO_W = 0.5
HEADER_H = 0.8           # inches for condition header
TITLE_H  = 0.7
CBAR_H   = 0.25
BOTTOM_PAD = 0.5
GAP      = 0.5           # gap between compartments

# Max proteins → figure height
max_n     = max(d['n'] for d in compartment_data)
heatmap_heights = {d['label']: d['n'] * ROW_H_IN for d in compartment_data}
max_hm_h  = max(heatmap_heights.values())
fig_h     = max_hm_h + HEADER_H + TITLE_H + BOTTOM_PAD + 1.0

n_rep_cols = len(compartment_data[0]['all_rep_cols'])
comp_w     = DENDRO_W + LABEL_W + n_rep_cols * COL_W
fig_w      = N_COMP * comp_w + (N_COMP - 1) * GAP + 1.0 + 0.6  # +cbar

fig = plt.figure(figsize=(fig_w, fig_h))

# Colorbar vabs across all datasets
all_vals = np.concatenate([d['mat'].ravel() for d in compartment_data])
vabs = float(np.nanpercentile(np.abs(all_vals[np.isfinite(all_vals)]), 96))
vabs = max(vabs, 0.5)

def fx(x): return x / fig_w
def fy(y): return y / fig_h

cbar_added = False

for ci, d in enumerate(compartment_data):
    label    = d['label']
    n        = d['n']
    mat      = d['mat']
    genes    = d['genes']
    groups   = d['groups']
    all_cols = d['all_rep_cols']
    Z        = d['linkage']
    hm_h     = heatmap_heights[label]
    n_cols   = mat.shape[1]

    # x positions for this compartment
    x0 = (ci * (comp_w + GAP) + 0.3) / fig_w

    dendro_left   = x0
    dendro_w_f    = DENDRO_W / fig_w
    label_left    = x0 + DENDRO_W / fig_w
    label_w_f     = LABEL_W / fig_w
    heat_left     = x0 + (DENDRO_W + LABEL_W) / fig_w
    heat_w_f      = (n_cols * COL_W) / fig_w

    # Top-align heatmaps
    heat_top      = (fig_h - TITLE_H - HEADER_H - 0.1) / fig_h
    heat_bottom_f = heat_top - hm_h / fig_h

    ax_heat   = fig.add_axes([heat_left, heat_bottom_f, heat_w_f, hm_h / fig_h])
    ax_labels = fig.add_axes([label_left, heat_bottom_f, label_w_f, hm_h / fig_h])
    ax_dendro = fig.add_axes([dendro_left, heat_bottom_f, dendro_w_f * 0.85, hm_h / fig_h])
    ax_header = fig.add_axes([heat_left, heat_top, heat_w_f, HEADER_H * 0.65 / fig_h])

    # Heatmap
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    im = ax_heat.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_heat.set_yticks([])
    ax_heat.set_xticks([])
    for sp in ax_heat.spines.values(): sp.set_visible(False)

    # White dividers between condition groups
    col_idx = 0
    for gname, gcol, gcols in groups:
        col_idx += len(gcols)
        if col_idx < n_cols:
            ax_heat.axvline(col_idx - 0.5, color='white', linewidth=2.0)

    # Dendrogram
    if Z is not None and n > 1:
        dendrogram(Z, ax=ax_dendro, orientation='left',
                   color_threshold=0, above_threshold_color='#777777',
                   no_labels=True, count_sort='ascending')
    ax_dendro.set_xlim(ax_dendro.get_xlim()[::-1])
    ax_dendro.axis('off')

    # Gene labels
    font_size = max(2.5, min(7, 350 / n))
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(-0.5, n - 0.5)
    ax_labels.invert_yaxis()
    for i, gene in enumerate(genes):
        ax_labels.text(0.95, i, gene, ha='right', va='center',
                       fontsize=font_size, fontfamily='Arial', color='#222222')
    ax_labels.axis('off')

    # Condition group header
    ax_header.set_xlim(0, n_cols)
    ax_header.set_ylim(0, 1)
    col_idx = 0
    tick_xs = []
    for gname, gcol, gcols in groups:
        nc = len(gcols)
        ax_header.add_patch(mpatches.FancyBboxPatch(
            (col_idx + 0.08, 0.1), nc - 0.16, 0.88,
            boxstyle='round,pad=0.02', facecolor=gcol, edgecolor='white', linewidth=1.5))
        ax_header.text(col_idx + nc / 2, 0.54, gname,
                       ha='center', va='center', fontsize=7.5,
                       fontweight='bold', fontfamily='Arial')
        for j, col in enumerate(gcols):
            short = str(col).split('_')[-1] if '_' in str(col) else str(col)
            ax_heat.text(col_idx + j, n + 0.5, short,
                         ha='center', va='top', fontsize=5.5, fontfamily='Arial',
                         transform=ax_heat.transData, rotation=55, color='#444444')
        col_idx += nc
    ax_header.axis('off')

    # Compartment title
    fig.text(heat_left + heat_w_f / 2,
             (fig_h - TITLE_H * 0.35) / fig_h,
             f'{label}   (n={n})',
             ha='center', va='center',
             fontsize=11, fontweight='bold', fontfamily='Arial')

    # Colorbar (only once, after last compartment)
    if ci == N_COMP - 1:
        cbar_left = heat_left + heat_w_f + 0.01
        cbar_w_f  = 0.012
        cbar_bot  = heat_bottom_f + hm_h * 0.3 / fig_h
        cbar_ht   = hm_h * 0.4 / fig_h
        ax_cbar   = fig.add_axes([cbar_left, cbar_bot, cbar_w_f, cbar_ht])
        cb = fig.colorbar(im, cax=ax_cbar)
        cb.set_label('Z-score', fontsize=8, fontfamily='Arial')
        cb.ax.tick_params(labelsize=7)

# Page title
fig.text(0.5, (fig_h - 0.18) / fig_h,
         'Subcellular proteome — union of significant proteins (corrected p > 3.3, |FC| > 0.5 in ≥1 condition)',
         ha='center', va='center', fontsize=10, fontfamily='Arial',
         color='#444444', style='italic')

fig.patch.set_facecolor('white')
plt.savefig(OUT, dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT}')
