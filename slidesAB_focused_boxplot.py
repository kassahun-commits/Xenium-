#!/usr/bin/env python3
"""
Box-and-whisker grid of the 47 focused-dot-plot genes.

One subplot per gene. Within each subplot:
  * x-axis: 7 condition groups
  * each group has Neuron + Astrocyte boxes side by side (dodged hue)
  * black dots overlay per-sample (per-mouse) mean expression of cells
    in that (gene, group, celltype) — so you can see the n=2-3 reps.
"""
from __future__ import annotations

import argparse
import logging
import sys
from datetime import date
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

LABEL_MAP = {
    'Chronic EtOH 1': ('ChronicEtOH', 1),
    'MCT1i+EtOH 1':   ('EtOH_MCT1i', 1),
    'MCT1i+Veh 1':    ('H2O_MCT1i',  1),
    'MAT2A CM 1':     ('MAT2A_CM',   1),
    'MAT2A CM 2':     ('MAT2A_CM',   2),
    'MAT2A CM 3':     ('MAT2A_CM',   3),
    'MAT2A OE 1':     ('MAT2A_OE',   1),
    'MAT2A OE 2':     ('MAT2A_OE',   2),
    'MAT2A OE 3':     ('MAT2A_OE',   3),
    'Veh+EtOH 1':     ('EtOH_veh',   1),
    'veh +EtOH 2':    ('EtOH_veh',   2),
    'veh+H20 1':      ('H2O_veh',    1),
    'EtOH +MCT1i 2':  ('EtOH_MCT1i', 2),
    'EtOH +MCT1i 3':  ('EtOH_MCT1i', 3),
    'H20+veh 2':      ('H2O_veh',    2),
    'H20 +MCT1i 2':   ('H2O_MCT1i',  2),
    'Chronic EtOH 2': ('ChronicEtOH', 2),
    'chronic EtOH 3': ('ChronicEtOH', 3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
          'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
CELLTYPES = ['Neuron', 'Astrocyte']
CT_COLORS = {'Neuron': '#1f77b4', 'Astrocyte': '#d62728'}


def load_polygons(csv):
    df = pd.read_csv(csv, comment='#')
    out = {}
    for name, sub in df.groupby('Selection', sort=False):
        out[name] = sub[['X', 'Y']].values.astype(float)
    return out


def assign(cells, polygons):
    out = pd.Series(np.full(len(cells), np.nan, dtype=object), index=cells.index)
    pts = cells[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        inside = MplPath(xy).contains_points(pts)
        out.values[inside & ~out.notna().values] = name
    return out


def process_slide(slide_dir, ann_csv, label):
    cells = pd.read_csv(slide_dir / 'cells.csv.gz').set_index('cell_id')
    polys = load_polygons(ann_csv)
    cells['roi'] = assign(cells, polys)
    cells['group'] = cells['roi'].map(lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)
    cells['slide'] = label
    cells_in = cells.dropna(subset=['group']).copy()
    adata = sc.read_10x_h5(str(slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    adata.obs_names = [f'{label}::{x}' for x in adata.obs_names]
    cells_in.index = [f'{label}::{x}' for x in cells_in.index]
    common = adata.obs_names.intersection(cells_in.index)
    adata = adata[common].copy()
    adata.obs = cells_in.loc[common, ['roi', 'group', 'replicate', 'sample_id', 'slide']].copy()
    return adata


def boxplot_panel(ax, gene_df, groups, celltypes, ct_colors, gene_name):
    """Draw one panel of boxes for one gene.

    Handles the zero-inflated case gracefully: if a gene's 90th percentile
    cell-level expression is < 0.5 (log-norm units), most cells are zeros
    and the box would collapse to an invisible line at y=0 — so we switch
    to a strip plot (colored dots = per-mouse means; no boxes) and label
    the panel as low-expression.
    """
    width = 0.32
    offsets = {'Neuron': -0.20, 'Astrocyte': 0.20}

    # Decide regime: standard box, or strip plot for sparse data
    all_vals = gene_df['expr'].values if len(gene_df) else np.array([0.0])
    q90 = float(np.percentile(all_vals, 90)) if len(all_vals) else 0.0
    use_strip = q90 < 0.5  # most cells are zero — boxes would be invisible

    y_max = 0.0
    for i, grp in enumerate(groups):
        for ct in celltypes:
            d = gene_df[(gene_df['group'] == grp) & (gene_df['celltype'] == ct)]
            if len(d) == 0:
                continue
            pos = i + offsets[ct]
            sample_means = d.groupby('sample_id')['expr'].mean()
            y_max = max(y_max, float(sample_means.max()) if len(sample_means) else 0.0)

            if not use_strip:
                ax.boxplot([d['expr'].values],
                           positions=[pos], widths=width,
                           patch_artist=True, showfliers=False,
                           medianprops=dict(color='black', lw=0.8),
                           whiskerprops=dict(color='black', lw=0.6),
                           capprops=dict(color='black', lw=0.6),
                           boxprops=dict(facecolor=ct_colors[ct],
                                         edgecolor='black', linewidth=0.5,
                                         alpha=0.6))
                # Per-mouse means as black dots on top of box
                jitter = np.random.uniform(-0.05, 0.05, size=len(sample_means))
                ax.scatter([pos] * len(sample_means) + jitter, sample_means.values,
                           c='black', s=10, zorder=5,
                           edgecolor='white', linewidth=0.4)
                y_max = max(y_max, float(np.percentile(d['expr'].values, 95)))
            else:
                # Strip plot: colored dots, no box
                jitter = np.random.uniform(-0.07, 0.07, size=len(sample_means))
                ax.scatter([pos] * len(sample_means) + jitter, sample_means.values,
                           c=ct_colors[ct], s=22, zorder=4,
                           edgecolor='black', linewidth=0.5, alpha=0.85)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha='right', fontsize=6)
    ax.set_title(gene_name, fontsize=9, fontweight='bold')
    ax.tick_params(axis='y', labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    # y-axis: always show 0; give a tiny pad below; expand to fit data
    ax.set_ylim(-0.05 * max(y_max, 0.05), max(y_max * 1.18, 0.1))
    if use_strip:
        ax.text(0.98, 0.97, 'low expr\n(per-mouse means)',
                transform=ax.transAxes, fontsize=5.5, va='top', ha='right',
                style='italic', color='#666',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--gene-source-csv', required=True, type=Path,
                    help='Focused dot plot source CSV (used to define gene list)')
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout)])

    # Genes: pull from the focused dot plot CSV in input order
    gene_df = pd.read_csv(args.gene_source_csv)
    genes = list(gene_df['gene'].drop_duplicates())
    logging.info('Genes to plot: %d', len(genes))

    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gs_present, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    score_cols = [f'score_{ct}' for ct in MARKERS]
    scores = adata.obs[score_cols]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_val = scores.max(axis=1)
    best_ct[best_val <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best_ct,
        categories=list(MARKERS.keys()) + ['Unclassified'])

    sub = adata[adata.obs['celltype'].astype(str).isin(CELLTYPES) &
                adata.obs['group'].astype(str).isin(GROUPS)].copy()
    logging.info('Cells in N/A + 7 groups: %s', f'{sub.shape[0]:,}')

    # Build the figure
    np.random.seed(0)
    n_genes = len(genes)
    n_cols = 7
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(2.4 * n_cols, 1.7 * n_rows),
                              squeeze=False)

    group_col = sub.obs['group'].astype(str).values
    ct_col = sub.obs['celltype'].astype(str).values
    sid_col = sub.obs['sample_id'].astype(str).values

    for i, gene in enumerate(genes):
        r, c = divmod(i, n_cols)
        ax = axes[r, c]
        if gene not in sub.var_names:
            ax.axis('off')
            continue
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        gdf = pd.DataFrame({'expr': v, 'group': group_col,
                            'celltype': ct_col, 'sample_id': sid_col})
        boxplot_panel(ax, gdf, GROUPS, CELLTYPES, CT_COLORS, gene)

    # Hide unused panels
    for j in range(n_genes, n_rows * n_cols):
        r, c = divmod(j, n_cols)
        axes[r, c].axis('off')

    # Legend (cell type)
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=CT_COLORS['Neuron'], edgecolor='black', alpha=0.6, label='Neuron'),
        Patch(facecolor=CT_COLORS['Astrocyte'], edgecolor='black', alpha=0.6, label='Astrocyte'),
    ]
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                      markerfacecolor='black', markersize=6,
                                      label='per-mouse mean'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=3, frameon=False, fontsize=10)

    fig.suptitle(
        f'{args.label}: focused gene panel — box plots across conditions in Neurons vs Astrocytes\n'
        f'(cell-level distributions; black dots = per-mouse means; n=2–3 per group)',
        fontsize=11, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.99])

    pdf_out = args.out_dir / f'{args.label}_FocusedBoxPlots_NeuronAstrocyte_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', pdf_out)


if __name__ == '__main__':
    main()
