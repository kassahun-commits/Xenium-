#!/usr/bin/env python3
"""
SlideB curated dot plot — alcohol + metabolism + activation gene panel.

One figure with two blocks: Neuron | Astrocyte. Each block has the 7
alcohol-relevant condition groups along its x-axis. Y-axis is ~100 genes
organized into categories (with banner labels on the right side and
separator lines between categories).

Reads cell type assignments from a previously saved per-cell CSV
(`SlideB_cells_with_roi.csv`) so we don't have to re-run the full
cell-typing pipeline.
"""
from __future__ import annotations

import argparse
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ----------------------------- Config -----------------------------

# Order on the x-axis (within each cell-type block)
GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
          'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']

CELLTYPES = ['Neuron', 'Astrocyte']

GENE_CATEGORIES = {
    'Cell-type identity': [
        'Rbfox3', 'Snap25', 'Syn1', 'Map2',
        'Gfap', 'Aqp4', 'Slc1a3', 'S100b', 'Aldh1l1',
    ],
    'Alcohol & acetaldehyde metab': [
        'Adh1', 'Adh5', 'Aldh2', 'Aldh1a1', 'Cyp2e1', 'Cat', 'Akr1a1',
    ],
    'Acetate / Ac-CoA metab': [
        'Acss1', 'Acss2', 'Acly', 'Acaca', 'Acacb',
    ],
    'Lactate / MCT axis': [
        'Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb',
    ],
    'Glucose transport': [
        'Slc2a1', 'Slc2a3', 'Slc2a4',
    ],
    'Glycolysis / TCA': [
        'Gapdh', 'Pgk1', 'Pkm', 'Hk1', 'Hk2', 'Pfkm', 'Aldoa',
        'Eno1', 'Eno2', 'Idh1', 'Idh2', 'Cs',
    ],
    'SAM / MAT2A axis': [
        'Mat1a', 'Mat2a', 'Mat2b', 'Ahcy', 'Mthfr', 'Cbs',
    ],
    'Glutathione / oxidative stress': [
        'Gclc', 'Gclm', 'Gpx1', 'Gpx4', 'Gss', 'Gstm1',
        'Sod1', 'Sod2', 'Hmox1', 'Nfe2l2',
    ],
    'Stress / UPR': [
        'Atf3', 'Atf4', 'Ddit3', 'Hspa1a', 'Hsp90aa1', 'Hspb1', 'Trp53',
    ],
    'Activity / IEGs': [
        'Fos', 'Jun', 'Junb', 'Arc', 'Egr1', 'Egr2',
        'Npas4', 'Homer1', 'Bdnf',
    ],
    'Glutamate signaling': [
        'Slc17a7', 'Slc17a6', 'Slc1a2',
        'Grin1', 'Grin2a', 'Grin2b',
        'Gria1', 'Gria2',
        'Grm1', 'Grm2', 'Grm5',
        'Camk2a', 'Camk2b',
    ],
    'GABA signaling': [
        'Gad1', 'Gad2', 'Gabra1', 'Gabra2', 'Gabbr1', 'Slc32a1',
    ],
    'Endocannabinoid': [
        'Cnr1', 'Faah', 'Mgll', 'Daglb',
    ],
    'Neuroinflammation': [
        'Tnf', 'Il1b', 'Il6', 'Il10', 'Tlr4', 'Nfkb1', 'Stat3',
    ],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--slide-dir', required=True, type=Path,
                   help='Xenium output folder with cell_feature_matrix.h5')
    p.add_argument('--cells-csv', required=True, type=Path,
                   help='Per-cell CSV with celltype/group columns (e.g. SlideB_cells_with_roi.csv)')
    p.add_argument('--out-dir', required=True, type=Path)
    p.add_argument('--label', default='SlideB')
    p.add_argument('--date', default=date.today().isoformat())
    return p.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load cell metadata (group, celltype) ----
    cells = pd.read_csv(args.cells_csv, index_col=0)
    needed = {'group', 'celltype'}
    if not needed.issubset(cells.columns):
        raise SystemExit(f'cells CSV missing columns {needed - set(cells.columns)}')
    print(f'Cell metadata rows: {len(cells):,}')

    # ---- Load counts ----
    print('Loading cell_feature_matrix.h5 ...')
    adata = sc.read_10x_h5(str(args.slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    print(f'AnnData: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes')

    common = adata.obs_names.intersection(cells.index)
    adata = adata[common].copy()
    adata.obs = cells.loc[common, ['group', 'celltype']].copy()

    # Subset to Neuron + Astrocyte; alcohol-relevant + MAT2A groups
    mask = (adata.obs['celltype'].astype(str).isin(CELLTYPES)) & \
           (adata.obs['group'].astype(str).isin(GROUPS))
    adata = adata[mask].copy()
    print(f'After subset (Neuron/Astrocyte; {len(GROUPS)} groups): '
          f'{adata.shape[0]:,} cells')

    # Normalize + log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Resolve gene availability
    all_listed = [g for cat in GENE_CATEGORIES.values() for g in cat]
    available = [g for g in all_listed if g in adata.var_names]
    missing = [g for g in all_listed if g not in adata.var_names]
    print(f'Genes listed: {len(all_listed)}; available in panel: {len(available)}; '
          f'missing: {len(missing)}')
    if missing:
        print(f'Missing: {missing}')

    # Build category -> available genes (in input order); drop empty categories
    cat_to_genes = {cat: [g for g in gs if g in available]
                    for cat, gs in GENE_CATEGORIES.items()}
    cat_to_genes = {k: v for k, v in cat_to_genes.items() if v}
    ordered_genes = [g for gs in cat_to_genes.values() for g in gs]
    print(f'Total genes to plot: {len(ordered_genes)}')

    # ---- Compute per-(group, celltype, gene) mean + fraction ----
    print('Computing per-group/celltype/gene statistics ...')
    group_col = adata.obs['group'].astype(str).values
    ct_col = adata.obs['celltype'].astype(str).values

    rows = []
    for gene in ordered_genes:
        X = adata[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        df = pd.DataFrame({'expr': v, 'group': group_col, 'celltype': ct_col})
        agg = (df.groupby(['group', 'celltype'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda x: float((x > 0).mean())),
                      n_cells=('expr', 'size'),
                      n_expr=('expr', lambda x: int((x > 0).sum())))
                 .reset_index())
        agg['gene'] = gene
        rows.append(agg)
    dot_df = pd.concat(rows, ignore_index=True)

    # Save source data
    csv_out = args.out_dir / f'{args.label}_CuratedDotPlot_NeuronAstrocyte_{args.date}.csv'
    dot_df.to_csv(csv_out, index=False)
    print(f'Wrote {csv_out}')

    # ---- Plot ----
    # x positions: 0-6 = Neuron block (7 groups); 8-14 = Astrocyte block (7 groups);
    #              gap at x=7 with a thin separator line.
    x_neuron = list(range(len(GROUPS)))
    x_astro  = [g + len(GROUPS) + 1 for g in range(len(GROUPS))]
    col_specs = [(g, 'Neuron', x_neuron[i]) for i, g in enumerate(GROUPS)] + \
                [(g, 'Astrocyte', x_astro[i]) for i, g in enumerate(GROUPS)]

    n_genes = len(ordered_genes)
    fig_w = 1.6 + (2 * len(GROUPS) + 1) * 0.45 + 1.8  # extra room on right for category labels
    fig_h = 1.5 + n_genes * 0.20
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmin = float(dot_df['mean_expr'].min())
    vmax = float(dot_df['mean_expr'].max())

    # Index lookups
    dot_idx = dot_df.set_index(['gene', 'group', 'celltype'])
    for i, gene in enumerate(ordered_genes):
        for grp, ct, x in col_specs:
            try:
                m = dot_idx.loc[(gene, grp, ct), 'mean_expr']
                f = dot_idx.loc[(gene, grp, ct), 'frac_expr']
            except KeyError:
                continue
            if pd.isna(m) or pd.isna(f):
                continue
            ax.scatter(x, i, s=(f * 220) + 5, c=[m], cmap='Reds',
                       vmin=vmin, vmax=vmax,
                       edgecolor='black', linewidth=0.4)

    # y-axis: gene names
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(ordered_genes, fontsize=7)
    ax.invert_yaxis()

    # x-axis: group labels twice (once per block)
    x_ticks = x_neuron + x_astro
    x_labels = GROUPS + GROUPS
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)

    # x-range
    ax.set_xlim(-0.8, x_astro[-1] + 0.8)
    ax.set_ylim(n_genes - 0.5, -2.0)

    # Block separator
    sep_x = len(GROUPS) - 0.5 + 1.0  # midway in the gap
    ax.axvline(sep_x, color='black', lw=0.8, alpha=0.5)

    # Block banner labels above
    ax.text((x_neuron[0] + x_neuron[-1]) / 2, -1.4, 'Neuron',
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text((x_astro[0] + x_astro[-1]) / 2, -1.4, 'Astrocyte',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Category separators (horizontal lines) + category banners on right
    label_x = x_astro[-1] + 1.0
    y_cur = 0
    cat_lines = []
    for cat, gs in cat_to_genes.items():
        n = len(gs)
        center = y_cur + (n - 1) / 2
        ax.text(label_x, center, cat, fontsize=8, va='center',
                ha='left', fontweight='bold')
        y_cur += n
        if y_cur < n_genes:
            cat_lines.append(y_cur - 0.5)
    for y in cat_lines:
        ax.axhline(y, color='black', lw=0.4, alpha=0.35)

    # Aesthetics
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap='Reds',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.18,
                        anchor=(1.0, 0.0))
    cbar.set_label('Mean log-norm expression', fontsize=8)

    # Size legend (manual, in top-right area outside axes)
    size_ax = fig.add_axes([0.92, 0.65, 0.08, 0.18])
    size_ax.axis('off')
    size_ax.set_title('% cells\nexpressing', fontsize=8, loc='left')
    for k, frac in enumerate([0.1, 0.25, 0.5, 0.75, 1.0]):
        size_ax.scatter(0.2, k, s=(frac * 220) + 5, c='lightgray',
                        edgecolor='black', linewidth=0.4)
        size_ax.text(0.55, k, f'{int(frac * 100)}%', va='center', fontsize=7)
    size_ax.set_xlim(0, 1); size_ax.set_ylim(-0.5, 4.5)

    fig.suptitle(
        f'{args.label}: curated alcohol / metabolism / activity gene panel '
        f'across conditions in Neurons vs Astrocytes',
        fontsize=11)
    pdf_out = args.out_dir / f'{args.label}_CuratedDotPlot_NeuronAstrocyte_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    print(f'Wrote {pdf_out}')


if __name__ == '__main__':
    main()
