import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmaps.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

COND_LABELS = ['Intoxication\nvs Naïve', 'Acute Withdrawal\nvs Naïve', 'Protracted Abstinence\nvs Naïve']

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_N-', 'Memb_{c}-', 'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-', 'Cyto_{c}-', 'all'),
    ('Chromatin',       'Chromatin',       'Chrom_N-','Chrom_{c}-','keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_N-',  'Nuc_{c}-',  'keep_review'),
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
    for sheet_label, excel_name, naive_pat, cond_tmpl, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
        if fmode == 'keep_review':
            df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        # Compute FC and significance for each condition
        fc_data    = {}
        corr_data  = {}
        for code in COND_CODES:
            naive_cols = [c for c in df.columns if naive_pat              in str(c)]
            cond_cols  = [c for c in df.columns if cond_tmpl.replace('{c}', code) in str(c)]
            fc, corr   = calc_stats(df, naive_cols, cond_cols)
            fc_data[code]   = fc
            corr_data[code] = corr

        # Union: significant in at least one condition
        sig_mask = pd.Series(False, index=df.index)
        for code in COND_CODES:
            fc   = fc_data[code]
            corr = corr_data[code]
            sig  = corr.notna() & fc.notna() & (corr > CORR_THRESH) & (fc.abs() > FC_THRESH)
            sig_mask = sig_mask | sig

        df_sig = df[sig_mask].reset_index(drop=True)
        n = len(df_sig)
        print(f'{sheet_label}: {n} proteins in union')

        if n == 0:
            print(f'  Skipping — no significant proteins')
            continue

        # Build FC matrix (rows=proteins, cols=conditions)
        fc_matrix = pd.DataFrame(index=df_sig.index)
        for code in COND_CODES:
            fc_matrix[code] = fc_data[code][sig_mask].values

        # Gene labels
        genes = df_sig['Gene symbol'].astype(str).tolist()

        # Hierarchical clustering of rows
        fc_vals = fc_matrix.values
        # Fill NaN with 0 for clustering only
        fc_for_clust = np.nan_to_num(fc_vals, nan=0.0)
        if n > 1:
            Z = linkage(fc_for_clust, method='ward', metric='euclidean')
            order = leaves_list(Z)
        else:
            order = [0]

        fc_ordered    = fc_vals[order, :]
        genes_ordered = [genes[i] for i in order]

        # ── Figure sizing ──
        col_width   = 2.5          # inches per condition column
        label_width = max(1.0, min(3.5, 200 / n))   # gene label column width
        dendro_w    = 1.2
        cbar_w      = 0.5
        fig_w       = dendro_w + label_width + col_width * 3 + cbar_w + 1.0

        row_h       = max(0.04, min(0.35, 14.0 / n))
        heatmap_h   = n * row_h
        top_pad     = 2.5         # title + column headers
        bot_pad     = 0.8
        fig_h       = heatmap_h + top_pad + bot_pad

        fig = plt.figure(figsize=(fig_w, fig_h))

        # GridSpec: [dendrogram | gene labels | heatmap | colorbar]
        gs = gridspec.GridSpec(
            1, 4,
            figure=fig,
            left=dendro_w / fig_w,
            right=1 - (cbar_w + 0.4) / fig_w,
            top=1 - top_pad / fig_h,
            bottom=bot_pad / fig_h,
            wspace=0.03,
            width_ratios=[dendro_w, label_width, col_width * 3, cbar_w],
        )

        ax_dendro = fig.add_axes([
            0.01,
            bot_pad / fig_h,
            (dendro_w - 0.1) / fig_w,
            heatmap_h / fig_h,
        ])
        ax_heat   = fig.add_axes([
            (dendro_w + label_width) / fig_w,
            bot_pad / fig_h,
            (col_width * 3) / fig_w,
            heatmap_h / fig_h,
        ])
        ax_labels = fig.add_axes([
            dendro_w / fig_w,
            bot_pad / fig_h,
            label_width / fig_w,
            heatmap_h / fig_h,
        ])
        ax_cbar   = fig.add_axes([
            (dendro_w + label_width + col_width * 3 + 0.15) / fig_w,
            bot_pad / fig_h + heatmap_h / fig_h * 0.2,
            0.18 / fig_w,
            heatmap_h / fig_h * 0.6,
        ])

        # Dendrogram
        if n > 1:
            dendrogram(Z, ax=ax_dendro, orientation='left',
                       color_threshold=0, above_threshold_color='#888888',
                       no_labels=True, count_sort='ascending')
        ax_dendro.axis('off')

        # Heatmap
        # Cap FC at ±2.5 for colour scale
        fc_disp = np.clip(fc_ordered, -2.5, 2.5)
        norm    = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
        im      = ax_heat.imshow(fc_disp, aspect='auto', cmap='RdBu_r', norm=norm,
                                 interpolation='nearest')
        ax_heat.set_xticks(range(3))
        ax_heat.set_xticklabels(COND_LABELS, fontsize=max(7, min(13, fig_w * 1.2)),
                                fontweight='bold', ha='center')
        ax_heat.xaxis.set_tick_params(length=0)
        ax_heat.tick_params(left=False, right=False, bottom=False)
        ax_heat.set_yticks([])
        for spine in ax_heat.spines.values():
            spine.set_visible(False)

        # Add thin grid lines between columns
        for x in [0.5, 1.5]:
            ax_heat.axvline(x, color='white', linewidth=1.5)

        # Gene labels
        font_size = max(3, min(9, 500 / n))
        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(-0.5, n - 0.5)
        for i, gene in enumerate(genes_ordered):
            ax_labels.text(0.95, i, gene, ha='right', va='center',
                           fontsize=font_size, fontfamily='Arial')
        ax_labels.axis('off')
        ax_labels.invert_yaxis()

        # Colorbar
        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label('Fold Change', fontsize=max(7, min(11, fig_w)),
                       fontfamily='Arial')
        cbar.ax.tick_params(labelsize=max(6, min(9, fig_w * 0.9)))

        # Title
        fig.text(0.5, 1 - 0.3 / fig_h,
                 f'{sheet_label}  —  Union of significant proteins  (n={n})',
                 ha='center', va='top',
                 fontsize=max(10, min(16, fig_w * 1.3)),
                 fontweight='bold', fontfamily='Arial')
        fig.text(0.5, 1 - 0.9 / fig_h,
                 f'Significant: |FC| > {FC_THRESH} and corrected p > {CORR_THRESH} in at least one condition',
                 ha='center', va='top',
                 fontsize=max(7, min(11, fig_w)),
                 color='#555555', style='italic', fontfamily='Arial')

        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Page saved.')

print(f'\nSaved: {OUT}')
