#!/usr/bin/env python3
"""
Xenium_May2026 — Slide B preliminary alcohol-focused analysis.

Pipeline:
  1. Load cells.csv.gz (centroids) + cell_feature_matrix.h5 from a Xenium output folder.
  2. Load polygon annotations exported from Xenium Explorer (one polygon per punch).
  3. Point-in-polygon: assign each cell to a ROI; map ROI label -> (group, replicate).
  4. QC filter, normalize total counts to 1e4, log1p.
  5. Marker-based cell typing (Neuron / Astrocyte / Oligodendrocyte / Microglia).
  6. Outputs: per-cell table, per-sample celltype counts, neuron:astrocyte ratio
     figure, MCT1-axis dotplot, QC figure. All figures vector-PDF with editable text.

Paths are passed via CLI (MEWS Lab rule 4: no hardcoded paths).
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

# Editable text in vector outputs (MEWS Lab rule 4.3)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ---------------------------------------------------------------------------
# Configuration (no paths — those come via CLI)
# ---------------------------------------------------------------------------

# Translate Xenium Explorer labels -> (group, replicate)
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
}

MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}

MCT_GENES = ['Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb']

# Alcohol-axis genes (requested by user)
ALCOHOL_AXIS_GENES = [
    'Adh5',    # alcohol dehydrogenase 5 / GSNOR
    'Adh1',    # alcohol dehydrogenase 1
    'Cat',     # catalase
    'Cyp2e1',  # CYP2E1, alcohol metabolism
    'Aldh2',   # acetaldehyde -> acetate
    'Acss1',   # acetyl-CoA synthetase short-chain 1
    'Acss2',   # acetyl-CoA synthetase short-chain 2
    'Slc2a1',  # GLUT1 (astrocyte/endothelial glucose transporter)
    'Slc2a3',  # GLUT3 (neuronal glucose transporter)
]

GROUP_ORDER = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
               'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
GROUP_COLORS = {
    'H2O_veh':     '#7f7f7f',
    'H2O_MCT1i':   '#17becf',
    'EtOH_veh':    '#ff7f0e',
    'EtOH_MCT1i':  '#2ca02c',
    'ChronicEtOH': '#d62728',
    'MAT2A_CM':    '#9467bd',
    'MAT2A_OE':    '#8c564b',
}
CT_COLORS = {
    'Neuron':          '#1f77b4',
    'Astrocyte':       '#d62728',
    'Oligodendrocyte': '#9467bd',
    'Microglia':       '#2ca02c',
    'Unclassified':    '#cccccc',
}


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--slide-dir', required=True, type=Path,
                   help='Xenium output folder (contains cells.csv.gz and cell_feature_matrix.h5)')
    p.add_argument('--annotations', required=True, type=Path,
                   help='Polygon annotations CSV from Xenium Explorer')
    p.add_argument('--out-dir', required=True, type=Path,
                   help='Output directory for figures + source data')
    p.add_argument('--label', default='SlideB', help='Slide label prefix for outputs')
    p.add_argument('--min-counts', type=int, default=10,
                   help='Min total transcripts per cell (QC)')
    p.add_argument('--min-cells', type=int, default=5,
                   help='Min cells expressing a gene (QC)')
    p.add_argument('--date', default=date.today().isoformat(),
                   help='Date tag for output filenames')
    return p.parse_args()


def setup_logging(out_dir: Path, label: str) -> Path:
    log_path = out_dir / f'{label}_log.txt'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[logging.FileHandler(log_path, mode='w'),
                  logging.StreamHandler(sys.stdout)],
    )
    return log_path


def load_polygons(csv_path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(csv_path, comment='#')
    if not {'Selection', 'X', 'Y'}.issubset(df.columns):
        raise ValueError(f'Annotations CSV missing required columns; got {df.columns.tolist()}')
    polygons = {}
    for name, sub in df.groupby('Selection', sort=False):
        polygons[name] = sub[['X', 'Y']].values.astype(float)
    return polygons


def aggregate_dotplot(adata, genes, row_col, col_col, row_order, col_order):
    """Compute per-(row,col) mean and fraction expressing for each gene.
    Returns long DataFrame: row, col, gene, mean_expr, frac_expr, n."""
    genes_present = [g for g in genes if g in adata.var_names]
    rows = []
    for gene in genes_present:
        X = adata[:, gene].X
        col_data = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        df = pd.DataFrame({'expr': col_data,
                           'row': adata.obs[row_col].astype(str).values,
                           'col': adata.obs[col_col].astype(str).values})
        agg = (df.groupby(['row', 'col'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda x: float((x > 0).mean())),
                      n=('expr', 'size'))
                 .reset_index())
        agg['gene'] = gene
        rows.append(agg)
    return pd.concat(rows, ignore_index=True), genes_present


def plot_dotplot(dot_df, genes, panels, row_order, title, out_pdf,
                 row_axis_label='', cbar_label='Mean log-norm expr'):
    """Faceted dot plot: one panel per value in `panels` (list).
    Each panel: rows = `row_order`, cols = `genes`. Dot size = frac_expr,
    color = mean_expr (shared scale across panels)."""
    import matplotlib.pyplot as plt
    vmin = float(dot_df['mean_expr'].min())
    vmax = float(dot_df['mean_expr'].max())
    n_panels = len(panels); n_g = len(genes); n_r = len(row_order)
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(1.6 + n_panels * (1.0 + n_g * 0.55),
                                      1.4 + n_r * 0.32),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]
    for ax, panel in zip(axes, panels):
        d = dot_df[dot_df['col'] == panel]
        mean_mat = (d.pivot(index='row', columns='gene', values='mean_expr')
                     .reindex(index=row_order, columns=genes))
        frac_mat = (d.pivot(index='row', columns='gene', values='frac_expr')
                     .reindex(index=row_order, columns=genes))
        for i, ri in enumerate(row_order):
            for j, gj in enumerate(genes):
                m = mean_mat.iat[i, j]; f = frac_mat.iat[i, j]
                if pd.isna(m) or pd.isna(f): continue
                ax.scatter(j, i, s=(f * 250) + 5, c=[m],
                           cmap='Reds', vmin=vmin, vmax=vmax,
                           edgecolor='black', linewidth=0.4)
        ax.set_xticks(range(n_g))
        ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(n_r))
        ax.set_yticklabels(row_order, fontsize=9)
        ax.set_title(panel, fontsize=10)
        ax.invert_yaxis()
        ax.set_xlim(-0.6, n_g - 0.4)
        ax.set_ylim(n_r - 0.4, -0.6)
        if row_axis_label and ax is axes[0]:
            ax.set_ylabel(row_axis_label, fontsize=9)
    sm = plt.cm.ScalarMappable(cmap='Reds',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_label(cbar_label, fontsize=8)
    fig.suptitle(title, fontsize=10)
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()


def assign_cells_to_rois(cells_df: pd.DataFrame,
                         polygons: dict[str, np.ndarray]) -> pd.Series:
    """Point-in-polygon. Returns ROI label or NaN per cell."""
    out = pd.Series(np.full(len(cells_df), np.nan, dtype=object),
                    index=cells_df.index, name='roi')
    pts = cells_df[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        path = MplPath(xy)
        inside = path.contains_points(pts)
        overlap = inside & out.notna().values
        if overlap.any():
            logging.warning(f'  Overlap: {overlap.sum()} cells already assigned, '
                            f'also inside "{name}" — keeping first assignment')
            inside = inside & ~out.notna().values
        out.values[inside] = name
    return out


# ---------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = setup_logging(out_dir, args.label)

    logging.info('=== %s alcohol prelim ===', args.label)
    logging.info('slide_dir:    %s', args.slide_dir)
    logging.info('annotations:  %s', args.annotations)
    logging.info('out_dir:      %s', out_dir)
    logging.info('QC: min_counts=%d, min_cells=%d', args.min_counts, args.min_cells)

    # ---- 1) Cells with centroids ----
    cells = pd.read_csv(args.slide_dir / 'cells.csv.gz').set_index('cell_id')
    logging.info('Cells loaded: %s', f'{len(cells):,}')

    # ---- 2) Polygons + ROI assignment ----
    polygons = load_polygons(args.annotations)
    logging.info('Polygons loaded: %d', len(polygons))
    cells['roi'] = assign_cells_to_rois(cells, polygons)
    n_assigned = cells['roi'].notna().sum()
    logging.info('Cells inside a ROI: %s / %s (%.1f%%)',
                 f'{n_assigned:,}', f'{len(cells):,}',
                 100.0 * n_assigned / len(cells))
    for roi, n in cells['roi'].value_counts(dropna=False).items():
        logging.info('  ROI=%s : %s', roi, f'{n:,}')

    # Label translation
    cells['group'] = cells['roi'].map(lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)

    # Warn about unmapped ROIs
    unmapped = set(cells['roi'].dropna().unique()) - set(LABEL_MAP.keys())
    if unmapped:
        logging.warning('Unmapped ROIs (will be excluded): %s', unmapped)

    cells_in = cells.dropna(subset=['group']).copy()
    logging.info('Cells assigned to a known group: %s', f'{len(cells_in):,}')

    # ---- 3) Load count matrix ----
    logging.info('Loading cell_feature_matrix.h5 ...')
    adata = sc.read_10x_h5(str(args.slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    logging.info('AnnData: %s cells x %s genes', f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

    common = adata.obs_names.intersection(cells_in.index)
    logging.info('Cells in both count matrix and ROIs: %s', f'{len(common):,}')
    adata = adata[common].copy()
    adata.obs = cells_in.loc[common,
        ['roi', 'group', 'replicate', 'sample_id', 'x_centroid', 'y_centroid']].copy()
    adata.obs['group'] = pd.Categorical(adata.obs['group'], categories=GROUP_ORDER, ordered=True)

    # ---- 4) QC ----
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    n_before = adata.shape[0]; g_before = adata.shape[1]
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    logging.info('After QC: %s cells (dropped %s) x %s genes (dropped %s)',
                 f'{adata.shape[0]:,}', f'{n_before - adata.shape[0]:,}',
                 f'{adata.shape[1]:,}', f'{g_before - adata.shape[1]:,}')

    # ---- 5) Normalize + log ----
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ---- 6) Marker scoring + cell type assignment ----
    score_cols = []
    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        if not gs_present:
            logging.warning('No markers found for %s', ct); continue
        sc.tl.score_genes(adata, gene_list=gs_present,
                          score_name=f'score_{ct}', random_state=0, n_bins=25)
        score_cols.append(f'score_{ct}')
        logging.info('  scored %s: used %d/%d markers: %s',
                     ct, len(gs_present), len(gs), gs_present)

    scores = adata.obs[score_cols]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_val = scores.max(axis=1)
    best_ct[best_val <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best_ct,
        categories=list(MARKERS.keys()) + ['Unclassified'])

    ct_counts = (adata.obs.groupby(['group', 'replicate', 'celltype'],
                                   observed=False)
                  .size().unstack(fill_value=0))
    logging.info('\nCell type counts per sample:\n%s', ct_counts.to_string())

    # ---- 7) Outputs ----
    L = args.label; D = args.date

    # 7a) per-cell table
    cell_out = out_dir / f'{L}_cells_with_roi.csv'
    adata.obs.to_csv(cell_out)
    logging.info('Wrote %s', cell_out)

    # 7b) celltype counts per sample
    counts_out = out_dir / f'{L}_celltype_counts_per_sample.csv'
    ct_counts.to_csv(counts_out)
    logging.info('Wrote %s', counts_out)

    # 7c) Neuron : Astrocyte ratio figure (all groups)
    df_n = ct_counts.copy()
    for c in MARKERS:
        if c not in df_n.columns: df_n[c] = 0
    df_n['n_total'] = df_n[list(MARKERS.keys()) + (['Unclassified'] if 'Unclassified' in df_n.columns else [])].sum(axis=1)
    df_n['frac_neuron']    = df_n['Neuron']    / df_n['n_total'].clip(lower=1)
    df_n['frac_astrocyte'] = df_n['Astrocyte'] / df_n['n_total'].clip(lower=1)
    df_n['neuron_per_astrocyte'] = df_n['Neuron'] / df_n['Astrocyte'].clip(lower=1)
    df_n_flat = df_n.reset_index()
    # Drop (group, replicate) combinations that have no cells (categorical produces all combos)
    df_n_flat = df_n_flat[df_n_flat['n_total'] > 0].copy()
    df_n_flat['sample_id'] = df_n_flat.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}", axis=1)
    df_n_flat = df_n_flat.sort_values(['group', 'replicate'])

    ratio_csv = out_dir / f'{L}_NeuronAstrocyteRatio_AllGroups_{D}.csv'
    df_n_flat.to_csv(ratio_csv, index=False)
    logging.info('Wrote %s', ratio_csv)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    samples = df_n_flat['sample_id'].tolist()
    bar_colors = [GROUP_COLORS.get(g, '#444') for g in df_n_flat['group']]

    # axes[0]: stacked fractions
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    ct_order = ['Neuron', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Unclassified']
    for ct in ct_order:
        if ct not in df_n_flat.columns: continue
        vals = (df_n_flat[ct] / df_n_flat['n_total']).values
        axes[0].bar(x, vals, bottom=bottom, width=0.8,
                    color=CT_COLORS.get(ct, '#888'), label=ct, edgecolor='white', linewidth=0.4)
        bottom += vals
    axes[0].set_xticks(x); axes[0].set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    axes[0].set_ylabel('Fraction of cells')
    axes[0].set_title(f'{L}: Cell-type composition per sample')
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)

    # axes[1]: neuron / astrocyte ratio
    axes[1].bar(x, df_n_flat['neuron_per_astrocyte'].values,
                color=bar_colors, edgecolor='black', linewidth=0.5)
    axes[1].set_xticks(x); axes[1].set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    axes[1].set_ylabel('Neuron / Astrocyte ratio')
    axes[1].set_title(f'{L}: Neuron : Astrocyte ratio')
    axes[1].axhline(1.0, ls='--', color='black', lw=0.6, alpha=0.5)

    # Build legend for group colors
    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor='black', label=g)
                      for g in GROUP_ORDER if g in df_n_flat['group'].astype(str).tolist()]
    axes[1].legend(handles=legend_handles, loc='center left',
                   bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)

    plt.tight_layout()
    ratio_pdf = out_dir / f'{L}_NeuronAstrocyteRatio_AllGroups_{D}.pdf'
    plt.savefig(ratio_pdf, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', ratio_pdf)

    # 7d) MCT1 axis dot plot — Neuron + Astrocyte, all groups
    mct_present = [g for g in MCT_GENES if g in adata.var_names]
    logging.info('MCT-axis genes present: %s', mct_present)

    main_cts = ['Neuron', 'Astrocyte']
    sub = adata[adata.obs['celltype'].isin(main_cts)].copy()
    dot_rows = []
    for gene in mct_present:
        col = np.asarray(sub[:, gene].X.todense()).ravel() if hasattr(sub[:, gene].X, 'todense') \
              else np.asarray(sub[:, gene].X).ravel()
        df_g = pd.DataFrame({'expr': col,
                             'group': sub.obs['group'].astype(str).values,
                             'ct': sub.obs['celltype'].astype(str).values})
        agg = (df_g.groupby(['group', 'ct'])
                   .agg(mean_expr=('expr', 'mean'),
                        frac_expr=('expr', lambda x: (x > 0).mean()),
                        n=('expr', 'size'))
                   .reset_index())
        agg['gene'] = gene
        dot_rows.append(agg)
    dot_df = pd.concat(dot_rows, ignore_index=True)
    mct_csv = out_dir / f'{L}_MCT1axis_dotplot_AllGroups_{D}.csv'
    dot_df.to_csv(mct_csv, index=False)
    logging.info('Wrote %s', mct_csv)

    # Dot-plot rendering: rows = group, cols = gene, one panel per cell type
    fig, axes = plt.subplots(1, len(main_cts),
                             figsize=(2.5 + len(mct_present) * 0.9,
                                      1.5 + len(GROUP_ORDER) * 0.32),
                             sharey=True)
    if len(main_cts) == 1:
        axes = [axes]

    vmin = dot_df['mean_expr'].min()
    vmax = dot_df['mean_expr'].max()
    for ax, ct in zip(axes, main_cts):
        d = dot_df[dot_df['ct'] == ct]
        mean_mat = (d.pivot(index='group', columns='gene', values='mean_expr')
                     .reindex(index=GROUP_ORDER, columns=mct_present))
        frac_mat = (d.pivot(index='group', columns='gene', values='frac_expr')
                     .reindex(index=GROUP_ORDER, columns=mct_present))
        for i, gi in enumerate(GROUP_ORDER):
            for j, gj in enumerate(mct_present):
                m = mean_mat.iat[i, j]; f = frac_mat.iat[i, j]
                if pd.isna(m) or pd.isna(f): continue
                ax.scatter(j, i, s=(f * 250) + 5, c=[m],
                           cmap='Reds', vmin=vmin, vmax=vmax,
                           edgecolor='black', linewidth=0.4)
        ax.set_xticks(range(len(mct_present)))
        ax.set_xticklabels(mct_present, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(GROUP_ORDER)))
        ax.set_yticklabels(GROUP_ORDER, fontsize=9)
        ax.set_title(ct, fontsize=10)
        ax.invert_yaxis()
        ax.set_xlim(-0.6, len(mct_present) - 0.4)
        ax.set_ylim(len(GROUP_ORDER) - 0.4, -0.6)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap='Reds',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_label('Mean log-norm expr', fontsize=8)
    fig.suptitle(f'{L}: MCT1-axis expression — dot size = fraction expressing',
                 fontsize=10)
    mct_pdf = out_dir / f'{L}_MCT1axis_dotplot_AllGroups_{D}.pdf'
    plt.savefig(mct_pdf, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', mct_pdf)

    # 7d.2) Alcohol-axis dot plot — Neuron + Astrocyte across groups
    sub_na = adata[adata.obs['celltype'].isin(main_cts)].copy()
    sub_na.obs['group_str'] = sub_na.obs['group'].astype(str)
    sub_na.obs['ct_str'] = sub_na.obs['celltype'].astype(str)
    etoh_df, etoh_present = aggregate_dotplot(
        sub_na, ALCOHOL_AXIS_GENES, row_col='group_str', col_col='ct_str',
        row_order=GROUP_ORDER, col_order=main_cts)
    etoh_df.rename(columns={'row': 'group', 'col': 'ct'}, inplace=True)
    etoh_csv = out_dir / f'{L}_AlcoholAxis_dotplot_AllGroups_{D}.csv'
    etoh_df.to_csv(etoh_csv, index=False)
    logging.info('Wrote %s', etoh_csv)
    logging.info('Alcohol-axis genes present: %s', etoh_present)
    # Plot via helper
    etoh_df_plot = etoh_df.rename(columns={'group': 'row', 'ct': 'col'})
    plot_dotplot(
        etoh_df_plot, genes=etoh_present, panels=main_cts,
        row_order=GROUP_ORDER,
        title=f'{L}: Alcohol-metabolism axis — dot size = fraction expressing',
        out_pdf=out_dir / f'{L}_AlcoholAxis_dotplot_AllGroups_{D}.pdf',
        row_axis_label='Group')
    logging.info('Wrote %s', out_dir / f'{L}_AlcoholAxis_dotplot_AllGroups_{D}.pdf')

    # 7d.3) Cell-type marker validation dot plot
    # Rows = assigned celltype; cols = marker gene; panel = marker lineage
    all_marker_genes = []
    gene_to_lineage = {}
    for ct, gs in MARKERS.items():
        present = [g for g in gs if g in adata.var_names]
        for g in present:
            if g not in gene_to_lineage:
                all_marker_genes.append(g)
                gene_to_lineage[g] = ct

    # Compute per (assigned_celltype, gene) mean and fraction
    val_rows = []
    cell_types_observed = [ct for ct in list(MARKERS.keys()) + ['Unclassified']
                            if (adata.obs['celltype'] == ct).any()]
    for gene in all_marker_genes:
        X = adata[:, gene].X
        col_data = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        for ct in cell_types_observed:
            mask = (adata.obs['celltype'].astype(str) == ct).values
            if mask.sum() == 0: continue
            val_rows.append({
                'celltype': ct, 'gene': gene,
                'lineage_of_marker': gene_to_lineage[gene],
                'mean_expr': float(col_data[mask].mean()),
                'frac_expr': float((col_data[mask] > 0).mean()),
                'n_cells': int(mask.sum()),
            })
    val_df = pd.DataFrame(val_rows)
    val_csv = out_dir / f'{L}_MarkerValidation_{D}.csv'
    val_df.to_csv(val_csv, index=False)
    logging.info('Wrote %s', val_csv)

    # Plot: single panel; rows = assigned celltype; cols = marker; group markers by lineage with separators
    fig, ax = plt.subplots(figsize=(2 + len(all_marker_genes) * 0.45,
                                    1 + len(cell_types_observed) * 0.4))
    vmin = val_df['mean_expr'].min(); vmax = val_df['mean_expr'].max()
    pivot_mean = val_df.pivot(index='celltype', columns='gene', values='mean_expr') \
                       .reindex(index=cell_types_observed, columns=all_marker_genes)
    pivot_frac = val_df.pivot(index='celltype', columns='gene', values='frac_expr') \
                       .reindex(index=cell_types_observed, columns=all_marker_genes)
    for i, ct in enumerate(cell_types_observed):
        for j, gene in enumerate(all_marker_genes):
            m = pivot_mean.iat[i, j]; f = pivot_frac.iat[i, j]
            if pd.isna(m) or pd.isna(f): continue
            ax.scatter(j, i, s=(f * 250) + 5, c=[m],
                       cmap='Reds', vmin=vmin, vmax=vmax,
                       edgecolor='black', linewidth=0.4)
    # Lineage separators on x
    cumulative = 0
    for ct, gs in MARKERS.items():
        present = [g for g in gs if g in adata.var_names]
        cumulative += len(present)
        if cumulative < len(all_marker_genes):
            ax.axvline(cumulative - 0.5, color='black', lw=0.5, alpha=0.4)
    # Lineage banner labels at top
    cumulative = 0
    for ct, gs in MARKERS.items():
        present = [g for g in gs if g in adata.var_names]
        if not present: continue
        center = cumulative + len(present) / 2 - 0.5
        ax.text(center, -0.9, ct, ha='center', va='bottom', fontsize=8, fontweight='bold')
        cumulative += len(present)
    ax.set_xticks(range(len(all_marker_genes)))
    ax.set_xticklabels(all_marker_genes, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(cell_types_observed)))
    ax.set_yticklabels(cell_types_observed, fontsize=9)
    ax.set_ylabel('Assigned cell type', fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(all_marker_genes) - 0.4)
    ax.set_ylim(len(cell_types_observed) - 0.4, -1.4)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('Mean log-norm expr', fontsize=8)
    ax.set_title(f'{L}: Marker validation — diagonal pattern = good cell typing\n'
                 f'(dot size = fraction expressing)', fontsize=10)
    plt.tight_layout()
    plt.savefig(out_dir / f'{L}_MarkerValidation_{D}.pdf', bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', out_dir / f'{L}_MarkerValidation_{D}.pdf')

    # 7f) Per-cell violin/strip plots of selected genes within Neuron / Astrocyte
    #     across the 5 alcohol-relevant groups
    EXPR_GENES = ['Adh5', 'Cat', 'Slc2a1', 'Slc2a3', 'Acss1', 'Acss2']
    ALC_GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i', 'ChronicEtOH']

    alc_mask = adata.obs['group'].astype(str).isin(ALC_GROUPS)
    adata_alc = adata[alc_mask].copy()

    for ct in ['Neuron', 'Astrocyte']:
        sub_ct = adata_alc[adata_alc.obs['celltype'].astype(str) == ct].copy()
        if sub_ct.shape[0] == 0:
            logging.warning('No %s cells in alcohol groups; skipping violin', ct)
            continue
        present = [g for g in EXPR_GENES if g in sub_ct.var_names]
        n_rows = (len(present) + 2) // 3
        fig, axes = plt.subplots(n_rows, 3, figsize=(13, 3.2 * n_rows),
                                  squeeze=False)
        # Build long-form dataframe
        long_rows = []
        for gene in present:
            X = sub_ct[:, gene].X
            v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
            df_g = pd.DataFrame({
                'expr': v,
                'group': sub_ct.obs['group'].astype(str).values,
                'sample_id': sub_ct.obs['sample_id'].astype(str).values,
            })
            df_g['gene'] = gene
            long_rows.append(df_g)
        long_df = pd.concat(long_rows, ignore_index=True)
        # Save source data
        viol_csv = out_dir / f'{L}_AlcoholGeneExpr_{ct}_{D}.csv'
        long_df.to_csv(viol_csv, index=False)
        logging.info('Wrote %s', viol_csv)
        # Plot
        for i, gene in enumerate(present):
            ax = axes.flat[i]
            d = long_df[long_df['gene'] == gene]
            # boxplot per group, ordered
            data_by_group = [d[d['group'] == g]['expr'].values for g in ALC_GROUPS]
            parts = ax.violinplot(data_by_group, showmeans=False,
                                  showmedians=True, widths=0.85)
            for pc, g in zip(parts['bodies'], ALC_GROUPS):
                pc.set_facecolor(GROUP_COLORS.get(g, '#888'))
                pc.set_edgecolor('black'); pc.set_alpha(0.7)
            # overlay sample-level means as black dots
            for j, g in enumerate(ALC_GROUPS, start=1):
                samples = d[d['group'] == g]['sample_id'].unique()
                for s in samples:
                    s_mean = d[(d['group'] == g) & (d['sample_id'] == s)]['expr'].mean()
                    ax.scatter(j + np.random.uniform(-0.05, 0.05), s_mean,
                               color='black', s=20, zorder=5, edgecolor='white', linewidth=0.5)
            ax.set_xticks(range(1, len(ALC_GROUPS) + 1))
            ax.set_xticklabels(ALC_GROUPS, rotation=30, ha='right', fontsize=8)
            ax.set_ylabel('log-norm expr', fontsize=8)
            ax.set_title(gene, fontsize=10)
        # Hide unused subplots
        for j in range(len(present), n_rows * 3):
            axes.flat[j].axis('off')
        fig.suptitle(f'{L}: {ct} — alcohol-axis gene expression across groups\n'
                     '(black dots = per-sample mean; n=1 per group for most)',
                     fontsize=11)
        plt.tight_layout()
        viol_pdf = out_dir / f'{L}_AlcoholGeneExpr_{ct}_{D}.pdf'
        plt.savefig(viol_pdf, bbox_inches='tight')
        plt.close()
        logging.info('Wrote %s', viol_pdf)

    # 7g) Exploratory DE within Neurons and Astrocytes
    #     EtOH_veh vs H2O_veh, and ChronicEtOH vs H2O_veh
    #     NOTE: cell-level Wilcoxon with n=1 per group is pseudoreplication —
    #     p-values are EXPLORATORY ONLY.
    DE_COMPARISONS = [
        ('EtOH_veh', 'H2O_veh'),
        ('ChronicEtOH', 'H2O_veh'),
        ('EtOH_MCT1i', 'EtOH_veh'),
    ]
    for ct in ['Neuron', 'Astrocyte']:
        sub_ct = adata_alc[adata_alc.obs['celltype'].astype(str) == ct].copy()
        for grp_test, grp_ref in DE_COMPARISONS:
            sub_de = sub_ct[sub_ct.obs['group'].astype(str).isin([grp_test, grp_ref])].copy()
            if sub_de.shape[0] < 50:
                logging.warning('Skipping DE %s vs %s in %s: too few cells (%d)',
                                grp_test, grp_ref, ct, sub_de.shape[0])
                continue
            sub_de.obs['de_group'] = sub_de.obs['group'].astype(str)
            try:
                sc.tl.rank_genes_groups(sub_de, 'de_group',
                                        method='wilcoxon',
                                        reference=grp_ref,
                                        n_genes=sub_de.shape[1])
            except Exception as e:
                logging.warning('DE %s vs %s in %s failed: %s',
                                grp_test, grp_ref, ct, e)
                continue
            result = sub_de.uns['rank_genes_groups']
            de_df = pd.DataFrame({
                'gene': [str(g) for g in result['names'][grp_test]],
                'logfc': result['logfoldchanges'][grp_test],
                'pval': result['pvals'][grp_test],
                'padj': result['pvals_adj'][grp_test],
                'score': result['scores'][grp_test],
            })
            de_df['comparison'] = f'{grp_test}_vs_{grp_ref}'
            de_df['celltype'] = ct
            de_df['n_cells_test'] = int((sub_de.obs['de_group'] == grp_test).sum())
            de_df['n_cells_ref'] = int((sub_de.obs['de_group'] == grp_ref).sum())
            de_csv = out_dir / f'{L}_DE_{ct}_{grp_test}_vs_{grp_ref}_{D}.csv'
            de_df.to_csv(de_csv, index=False)
            logging.info('Wrote %s (top10 up: %s)', de_csv,
                         de_df.head(10)['gene'].tolist())
            # Volcano plot
            fig, ax = plt.subplots(figsize=(6, 5))
            x = de_df['logfc'].values
            y = -np.log10(np.clip(de_df['padj'].values, 1e-300, 1))
            sig = (de_df['padj'] < 0.01) & (de_df['logfc'].abs() > 0.5)
            ax.scatter(x[~sig], y[~sig], s=6, c='lightgray', alpha=0.5,
                       edgecolor='none')
            ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=12, c='#d62728',
                       alpha=0.8, edgecolor='none', label=f'Up in {grp_test}')
            ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=12, c='#1f77b4',
                       alpha=0.8, edgecolor='none', label=f'Up in {grp_ref}')
            # Label top genes
            top_up = de_df[sig & (de_df['logfc'] > 0)].nlargest(8, 'logfc')
            top_dn = de_df[sig & (de_df['logfc'] < 0)].nsmallest(8, 'logfc')
            for _, r in pd.concat([top_up, top_dn]).iterrows():
                ax.annotate(r['gene'], (r['logfc'], -np.log10(max(r['padj'], 1e-300))),
                            fontsize=7, alpha=0.85)
            ax.axhline(-np.log10(0.01), ls='--', color='black', lw=0.5, alpha=0.4)
            ax.axvline(0.5, ls='--', color='black', lw=0.5, alpha=0.4)
            ax.axvline(-0.5, ls='--', color='black', lw=0.5, alpha=0.4)
            ax.set_xlabel(f'log2 fold-change ({grp_test} / {grp_ref})', fontsize=9)
            ax.set_ylabel('-log10 adj p-value', fontsize=9)
            ax.set_title(
                f'{L}: {ct} — {grp_test} vs {grp_ref}\n'
                f'EXPLORATORY (n=1 mouse per group; cell-level Wilcoxon)',
                fontsize=10)
            ax.legend(fontsize=7, frameon=False)
            plt.tight_layout()
            de_pdf = out_dir / f'{L}_DE_{ct}_{grp_test}_vs_{grp_ref}_{D}.pdf'
            plt.savefig(de_pdf, bbox_inches='tight')
            plt.close()
            logging.info('Wrote %s', de_pdf)

    # 7e) QC figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].hist(adata.obs['total_counts'], bins=80, color='#444')
    axes[0].set_xlabel('total_counts'); axes[0].set_ylabel('cells')
    axes[0].set_title('Counts per cell (post-QC)')
    axes[0].set_yscale('log')

    s_counts = adata.obs.groupby('sample_id', observed=True).size().sort_values(ascending=False)
    sc_colors = [GROUP_COLORS.get(s.rsplit('_', 1)[0], '#444') for s in s_counts.index]
    axes[1].bar(range(len(s_counts)), s_counts.values, color=sc_colors,
                edgecolor='black', linewidth=0.4)
    axes[1].set_xticks(range(len(s_counts)))
    axes[1].set_xticklabels(s_counts.index, rotation=45, ha='right', fontsize=8)
    axes[1].set_ylabel('cells'); axes[1].set_title('Cells per sample (post-QC)')

    ct_overall = adata.obs['celltype'].value_counts().reindex(
        list(MARKERS.keys()) + ['Unclassified']).fillna(0)
    axes[2].bar(range(len(ct_overall)), ct_overall.values,
                color=[CT_COLORS.get(c, '#888') for c in ct_overall.index],
                edgecolor='black', linewidth=0.4)
    axes[2].set_xticks(range(len(ct_overall)))
    axes[2].set_xticklabels(ct_overall.index, rotation=45, ha='right', fontsize=8)
    axes[2].set_ylabel('cells'); axes[2].set_title('Cell types overall')
    plt.tight_layout()
    qc_pdf = out_dir / f'{L}_QCmetrics_{D}.pdf'
    plt.savefig(qc_pdf, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', qc_pdf)

    logging.info('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
