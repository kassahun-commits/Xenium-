import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib.patches as mpatches

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmaps_v2.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# Pastel diverging colormap: soft blue → white → soft rose
CMAP = LinearSegmentedColormap.from_list(
    'pastel_div',
    ['#6BAED6', '#AED6F0', '#FFFFFF', '#F4A8A8', '#E06060'],
    N=256
)

# Condition group colours for header bands (pastel)
GROUP_COLORS = {
    'Naïve':                  '#D4E6F1',
    'Intoxication':           '#FAD7A0',
    'Acute withdrawal':       '#A9DFBF',
    'Protracted abstinence':  '#D7BDE2',
}

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_N-', 'Memb_I-', 'Memb_AW-', 'Memb_PA-', 'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-', 'Cyto_I-', 'Cyto_AW-', 'Cyto_PA-', 'all'),
    ('Chromatin',       'Chromatin',       'Chrom_N-','Chrom_I-','Chrom_AW-','Chrom_PA-','keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_N-',  'Nuc_I-',  'Nuc_AW-',  'Nuc_PA-',  'keep_review'),
]
COND_CODES = ['I', 'AW', 'PA']

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

with PdfPages(OUT) as pdf:
    for sheet_label, excel_name, pat_N, pat_I, pat_AW, pat_PA, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
        if fmode == 'keep_review':
            df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        # Identify replicate columns for each group
        groups = [
            ('Naïve',               [c for c in df.columns if pat_N  in str(c)]),
            ('Intoxication',        [c for c in df.columns if pat_I  in str(c)]),
            ('Acute withdrawal',    [c for c in df.columns if pat_AW in str(c)]),
            ('Protracted abstinence',[c for c in df.columns if pat_PA in str(c)]),
        ]
        all_rep_cols = [c for _, cols in groups for c in cols]

        # Find union of significant proteins
        sig_mask = pd.Series(False, index=df.index)
        for code, pat in [('I', pat_I), ('AW', pat_AW), ('PA', pat_PA)]:
            naive_cols = [c for c in df.columns if pat_N in str(c)]
            cond_cols  = [c for c in df.columns if pat  in str(c)]
            fc, corr   = calc_stats(df, naive_cols, cond_cols)
            sig = corr.notna() & fc.notna() & (corr > CORR_THRESH) & (fc.abs() > FC_THRESH)
            sig_mask = sig_mask | sig

        df_sig = df[sig_mask].reset_index(drop=True)
        n = len(df_sig)
        print(f'{sheet_label}: {n} proteins')

        if n == 0:
            continue

        # Raw replicate matrix, numeric
        mat = df_sig[all_rep_cols].apply(pd.to_numeric, errors='coerce').values.astype(float)

        # Z-score each row (per protein across all replicates)
        row_mean = np.nanmean(mat, axis=1, keepdims=True)
        row_std  = np.nanstd(mat,  axis=1, keepdims=True)
        row_std[row_std == 0] = 1
        mat_z = (mat - row_mean) / row_std

        # Hierarchical clustering
        mat_fill = np.where(np.isnan(mat_z), 0, mat_z)
        if n > 1:
            Z = linkage(mat_fill, method='ward')
            order = leaves_list(Z)
        else:
            order = [0]

        mat_ordered   = mat_z[order, :]
        genes_ordered = df_sig['Gene symbol'].astype(str).iloc[order].tolist()

        # ── Figure layout ──
        n_cols      = len(all_rep_cols)
        col_w       = 0.35           # inches per replicate column
        label_w     = max(0.8, min(3.0, 180 / n))
        dendro_w    = 0.8
        cbar_w      = 0.25
        gap         = 0.2
        fig_w       = dendro_w + label_w + n_cols * col_w + cbar_w + gap + 0.5

        row_h       = max(0.03, min(0.30, 12.0 / n))
        heatmap_h   = n * row_h
        header_h    = 1.2            # condition group header
        top_pad     = 2.2
        bot_pad     = 0.6
        fig_h       = heatmap_h + top_pad + bot_pad

        fig = plt.figure(figsize=(fig_w, fig_h))

        # Axes positions (as fractions)
        def fx(x): return x / fig_w
        def fy(y): return y / fig_h

        heat_left   = fx(dendro_w + label_w)
        heat_bottom = fy(bot_pad)
        heat_w      = fx(n_cols * col_w)
        heat_h      = fy(heatmap_h)

        ax_dendro = fig.add_axes([fx(0.02), fy(bot_pad), fx(dendro_w - 0.1), heat_h])
        ax_labels = fig.add_axes([fx(dendro_w), fy(bot_pad), fx(label_w), heat_h])
        ax_heat   = fig.add_axes([heat_left, heat_bottom, heat_w, heat_h])
        ax_header = fig.add_axes([heat_left, heat_bottom + heat_h, heat_w, fy(header_h * 0.55)])
        ax_cbar   = fig.add_axes([heat_left + heat_w + fx(gap),
                                   heat_bottom + heat_h * 0.25,
                                   fx(cbar_w), heat_h * 0.5])

        # Dendrogram
        if n > 1:
            from scipy.cluster.hierarchy import dendrogram as dend
            dend(Z, ax=ax_dendro, orientation='left', color_threshold=0,
                 above_threshold_color='#888888', no_labels=True, count_sort='ascending')
        ax_dendro.axis('off')

        # Heatmap
        vabs = min(2.5, np.nanpercentile(np.abs(mat_ordered), 98))
        norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        im   = ax_heat.imshow(mat_ordered, aspect='auto', cmap=CMAP, norm=norm,
                              interpolation='nearest')
        ax_heat.set_yticks([])
        ax_heat.set_xticks([])
        for sp in ax_heat.spines.values(): sp.set_visible(False)

        # White dividers between condition groups
        col_idx = 0
        for gname, gcols in groups:
            col_idx += len(gcols)
            if col_idx < n_cols:
                ax_heat.axvline(col_idx - 0.5, color='white', linewidth=2.5)

        # Gene labels
        font_size = max(3, min(8, 450 / n))
        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(-0.5, n - 0.5)
        ax_labels.invert_yaxis()
        for i, gene in enumerate(genes_ordered):
            ax_labels.text(0.95, i, gene, ha='right', va='center',
                           fontsize=font_size, fontfamily='Arial')
        ax_labels.axis('off')

        # Condition group header
        ax_header.set_xlim(0, n_cols)
        ax_header.set_ylim(0, 1)
        col_idx = 0
        for gname, gcols in groups:
            nc = len(gcols)
            color = GROUP_COLORS[gname]
            ax_header.add_patch(mpatches.FancyBboxPatch(
                (col_idx + 0.05, 0.05), nc - 0.1, 0.9,
                boxstyle='round,pad=0.02', facecolor=color, edgecolor='white', linewidth=1.5
            ))
            ax_header.text(col_idx + nc / 2, 0.52, gname,
                           ha='center', va='center',
                           fontsize=max(6, min(11, fig_w * 0.7)),
                           fontweight='bold', fontfamily='Arial')
            # Individual replicate tick labels (short)
            for j, col in enumerate(gcols):
                short = col.split('_')[-1] if '_' in col else col
                ax_heat.text(col_idx + j, -0.8, short,
                             ha='center', va='top',
                             fontsize=max(4, min(7, fig_w * 0.5)),
                             fontfamily='Arial',
                             transform=ax_heat.get_xaxis_transform(),
                             rotation=60)
            col_idx += nc
        ax_header.axis('off')

        # Colorbar
        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label('Z-score', fontsize=max(6, min(10, fig_w * 0.8)), fontfamily='Arial')
        cbar.ax.tick_params(labelsize=max(5, min(8, fig_w * 0.7)))

        # Title
        fig.text(0.5, 1 - 0.25 / fig_h,
                 f'{sheet_label}  —  n = {n} proteins (significant in ≥1 condition)',
                 ha='center', va='top',
                 fontsize=max(10, min(15, fig_w * 1.1)),
                 fontweight='bold', fontfamily='Arial')
        fig.text(0.5, 1 - 0.75 / fig_h,
                 'Rows z-scored across all replicates  |  '
                 f'|FC| > {FC_THRESH} and corrected p > {CORR_THRESH} in at least one condition',
                 ha='center', va='top',
                 fontsize=max(6, min(10, fig_w * 0.75)),
                 color='#555555', style='italic', fontfamily='Arial')

        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved.')

print(f'\nSaved: {OUT}')
