"""
Heatmap per compartment: Intoxication / Acute Withdrawal / Protracted Abstinence
  - Rows  : proteins significant (corrected p > 3.3, |FC| > 0.5) in >= 1 condition
  - Cols  : ΔIntox, ΔAW, ΔPA  (fold change vs Naïve, log2 LFQ space)
  - Color : raw FC values (red = up, blue = down), clipped at ±3
  - Filter: Soluble Nuclear → all rows (Keep+Review+Exclude)
            Chromatin       → Keep + Review only
            Membrane/Cytosol → no Filter column, use all
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmap_Conditions_AllCompartments.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
CLIP        = 3.0    # clip FC values at ±3 for colour scale

CMAP = LinearSegmentedColormap.from_list('rwb',
    ['#2166AC','#74ADD1','#ABD9E9','#FFFFFF','#FDDBC7','#F4A582','#D6604D','#B2182B'], N=256)

SHEETS = [
    # (display_label, excel_sheet, prefix, filter_mode)
    ('Membrane',        'Membrane',        'Memb_',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'all'),   # includes Exclude
]

COND_SUFFIXES = {
    'Naive':     ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Intox':     ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'AW':        ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
    'PA':        ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
CONDITIONS = ['Intox', 'AW', 'PA']

def get_cols(df, prefix, suffixes):
    return [c for c in df.columns if any(prefix + s in str(c) for s in suffixes)]

def calc_fc_and_corr(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pv    = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pv.append(stats.ttest_ind(n, c, equal_var=False)[1]
                  if len(n) >= 2 and len(c) >= 2 else np.nan)
    pv   = pd.Series(pv, index=df.index)
    valid = pv.notna()
    ranks = pv[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pv[valid] * valid.sum() / ranks))
    return fc, corr

with PdfPages(OUT) as pdf:

    for disp_label, sheet, prefix, fmode in SHEETS:
        print(f'Loading {disp_label}...', end=' ', flush=True)
        raw = pd.read_excel(FILE, sheet_name=sheet, header=0)

        # Filter rows
        if fmode == 'keep_review' and 'Filter' in raw.columns:
            df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)  # all rows incl. Exclude

        print(f'{len(df)} proteins')

        # Column groups
        naive_cols = get_cols(df, prefix, COND_SUFFIXES['Naive'])
        cond_cols  = {c: get_cols(df, prefix, COND_SUFFIXES[c]) for c in CONDITIONS}

        # FC and corrected p per condition
        fc_dict   = {}
        corr_dict = {}
        for cond in CONDITIONS:
            fc, corr = calc_fc_and_corr(df, naive_cols, cond_cols[cond])
            fc_dict[cond]   = fc
            corr_dict[cond] = corr

        # Significant in at least one condition
        sig_mask = pd.Series(False, index=df.index)
        for cond in CONDITIONS:
            sig_mask |= (corr_dict[cond] > CORR_THRESH) & (fc_dict[cond].abs() > FC_THRESH)

        sig_df = df.loc[sig_mask].copy()
        genes  = sig_df['Gene symbol'].astype(str).values
        n_sig  = len(sig_df)
        print(f'  → {n_sig} proteins significant in ≥1 condition')

        # Build FC matrix (n_sig × 3)
        mat = np.column_stack([
            fc_dict[c].loc[sig_mask].values for c in CONDITIONS
        ]).astype(float)

        # Remove rows with all-NaN
        valid_rows = ~np.all(np.isnan(mat), axis=1)
        mat   = mat[valid_rows]
        genes = genes[valid_rows]
        n_sig = len(mat)

        # Fill remaining NaN with 0 for clustering (no change = neutral)
        mat_fill = np.where(np.isnan(mat), 0, mat)

        # Hierarchical clustering (Ward on rows)
        if n_sig > 1:
            dist  = pdist(mat_fill, metric='euclidean')
            Z     = linkage(dist, method='ward')
            order = leaves_list(Z)
        else:
            order = np.arange(n_sig)

        mat_sorted   = mat[order]
        genes_sorted = genes[order]

        # Clip for colour scale
        mat_clipped = np.clip(mat_sorted, -CLIP, CLIP)

        # ── Plot ──────────────────────────────────────────────────────────────
        # Dynamic figure height based on n_sig
        row_h  = max(0.012, min(0.06, 400 / n_sig))   # inches per row
        fig_h  = max(6, min(28, n_sig * row_h + 3))
        fig_w  = 7

        fig = plt.figure(figsize=(fig_w, fig_h))
        fig.patch.set_facecolor('white')

        gs = gridspec.GridSpec(1, 2, figure=fig,
                               width_ratios=[1, 0.04],
                               wspace=0.04,
                               left=0.22, right=0.88,
                               top=0.93, bottom=0.06)

        ax     = fig.add_subplot(gs[0, 0])
        ax_cb  = fig.add_subplot(gs[0, 1])

        im = ax.imshow(mat_clipped, aspect='auto', cmap=CMAP,
                       vmin=-CLIP, vmax=CLIP, interpolation='nearest')

        # X-axis: condition labels
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['Intoxication', 'Acute\nWithdrawal', 'Protracted\nAbstinence'],
                           fontsize=11, fontweight='bold')
        ax.xaxis.set_tick_params(length=0)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')

        # Y-axis: gene labels only if few enough to read
        if n_sig <= 80:
            ax.set_yticks(range(n_sig))
            ax.set_yticklabels(genes_sorted, fontsize=max(4, 8 - n_sig // 20))
        else:
            ax.set_yticks([])
            ax.set_ylabel(f'{n_sig} proteins', fontsize=11, labelpad=8)

        # Vertical dividers between conditions
        for x in [0.5, 1.5]:
            ax.axvline(x, color='white', linewidth=2)

        ax.set_title(f'{disp_label}\n(n = {n_sig} proteins, sig. in ≥1 condition)',
                     fontsize=13, fontweight='bold', pad=14)

        # Colourbar
        cb = plt.colorbar(im, cax=ax_cb)
        cb.set_label('Fold Change\nvs Naïve', fontsize=10, labelpad=8)
        cb.ax.tick_params(labelsize=9)
        cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
        cb.set_ticklabels(['≤−3', '−2', '−1', '0', '+1', '+2', '≥+3'])

        pdf.savefig(fig, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved page for {disp_label}')

print(f'\nDone → {OUT}')
