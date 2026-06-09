#!/usr/bin/env python3
"""
Run the 10x Genomics "Analyze a single Xenium sample using Python" workshop
pipeline on the WHOLE-CELL region-V1 data (Slides A + B), starting from the V1
ROI coordinates we already have.

Faithful to the 10x workshop notebook
(github.com/10XGenomics/analysis_guides/blob/main/
 Xenium_workshop_single_sample_analysis.ipynb):

  QC        filter_cells(min_counts=20); filter_cells(max_counts=q98);
            filter_genes(min_cells=100)
  HVG       highly_variable_genes(flavor="seurat_v3", n_top_genes=2000)
  Norm      normalize_total (median) -> log1p -> layers["lognorm"]
            -> scale(zero_center=False, max_value=10)
  Cluster   pca(n_comps=30) -> neighbors(metric="cosine")
            -> leiden(resolution=0.5, flavor="igraph") -> umap
  Markers   rank_genes_groups(groupby="leiden", layer="lognorm") + dotplot
  Spatial   sq.gr.spatial_neighbors(delaunay=True) [PER SLIDE via library_key]
            -> sq.gr.nhood_enrichment(cluster_key="leiden") -> heatmap

Differences from the notebook (documented, intentional):
  - We build the AnnData straight from the raw output files and restrict to the
    V1 ROIs with our own point-in-polygon assignment (no spatialdata/zarr).
  - Combined V1 = both slides. The workshop is single-sample and does NO batch
    correction, so leiden/UMAP here can carry slide batch effects (caveat). The
    spatial graph is built PER SLIDE (library_key="slide") so neighbours are
    never linked across physical sections.

No hardcoded paths (CLI). Editable-text PDFs (fonttype 42); source CSVs saved.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import squidpy as sq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


def seurat_v3_hvg(adata, n_top_genes, layer='counts'):
    """Faithful reproduction of scanpy's flavor='seurat_v3' HVG selection.

    Uses statsmodels lowess in place of skmisc.loess: scikit-misc ships no
    numpy-2-compatible wheel for Python 3.9 and there is no Fortran compiler on
    this machine to build it from source. The mean-variance fit (loess span 0.3
    -> lowess frac 0.3) and the clipped standardized-variance ranking otherwise
    match scanpy._highly_variable_genes_seurat_v3 exactly (single batch).
    Sets adata.var['highly_variable'] + means/variances/variances_norm.
    """
    from scipy import sparse
    from statsmodels.nonparametric.smoothers_lowess import lowess

    X = adata.layers[layer] if layer in adata.layers else adata.X
    X = sparse.csc_matrix(X).astype(np.float64)
    N = X.shape[0]
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = (mean_sq - mean ** 2) * (N / (N - 1))           # ddof=1
    not_const = var > 0

    estimat_var = np.zeros(X.shape[1], dtype=np.float64)
    y = np.log10(var[not_const])
    x = np.log10(mean[not_const])
    estimat_var[not_const] = lowess(y, x, frac=0.3, it=0, return_sorted=False)
    reg_std = np.sqrt(10 ** estimat_var)

    # Seurat vst: clip raw counts at reg_std*sqrt(N)+mean, take variance of the
    # clipped, standardized values. Zeros clip to 0 -> contribute 0 to the sums.
    clip_val = reg_std * np.sqrt(N) + mean
    indptr, data = X.indptr, X.data
    sum_c = np.zeros(X.shape[1])
    sum_c2 = np.zeros(X.shape[1])
    for j in range(X.shape[1]):
        col = np.minimum(data[indptr[j]:indptr[j + 1]], clip_val[j])
        sum_c[j] = col.sum()
        sum_c2[j] = (col * col).sum()
    with np.errstate(invalid='ignore', divide='ignore'):
        norm_var = (sum_c2 - 2 * mean * sum_c + N * mean ** 2) / \
                   ((N - 1) * reg_std ** 2)
    norm_var[~not_const] = 0.0

    order = np.argsort(norm_var)[::-1]
    hv = np.zeros(X.shape[1], dtype=bool)
    hv[order[:n_top_genes]] = True
    adata.var['highly_variable'] = hv
    adata.var['means'] = mean
    adata.var['variances'] = var
    adata.var['variances_norm'] = norm_var
    return adata

# ---- V1 ROI label map (verbatim from build_PI_master_NvsA_lfc0.5.py) ----
LABEL_MAP = {
    'Chronic EtOH 1': ('ChronicEtOH', 1), 'MCT1i+EtOH 1': ('EtOH_MCT1i', 1),
    'MCT1i+Veh 1': ('H2O_MCT1i', 1), 'MAT2A CM 1': ('MAT2A_CM', 1),
    'MAT2A CM 2': ('MAT2A_CM', 2), 'MAT2A CM 3': ('MAT2A_CM', 3),
    'MAT2A OE 1': ('MAT2A_OE', 1), 'MAT2A OE 2': ('MAT2A_OE', 2),
    'MAT2A OE 3': ('MAT2A_OE', 3), 'Veh+EtOH 1': ('EtOH_veh', 1),
    'veh +EtOH 2': ('EtOH_veh', 2), 'veh+H20 1': ('H2O_veh', 1),
    'EtOH +MCT1i 2': ('EtOH_MCT1i', 2), 'EtOH +MCT1i 3': ('EtOH_MCT1i', 3),
    'H20+veh 2': ('H2O_veh', 2), 'H20 +MCT1i 2': ('H2O_MCT1i', 2),
    'Chronic EtOH 2': ('ChronicEtOH', 2), 'chronic EtOH 3': ('ChronicEtOH', 3),
    'EtOH+MCT1i 1': ('EtOH_MCT1i', 1), 'H20+MCT1i 1': ('H2O_MCT1i', 1),
    'H20 +MCT1i 1': ('H2O_MCT1i', 1), 'H20 +veh 1': ('H2O_veh', 1),
    'EtOH +veh 1': ('EtOH_veh', 1), 'APP 1': ('APP', 1), 'APP 2': ('APP', 2),
    'APP 3': ('APP', 3),
}

# 10x workshop QC parameters
MIN_COUNTS = 20
MAX_COUNTS_Q = 0.98
MIN_CELLS = 100
N_HVG = 2000
N_PCS = 30
LEIDEN_RES = 0.5


def load_polygons(csv):
    df = pd.read_csv(csv, comment='#')
    return {name: sub[['X', 'Y']].values.astype(float)
            for name, sub in df.groupby('Selection', sort=False)}


def assign(cells, polygons):
    out = pd.Series(np.full(len(cells), np.nan, dtype=object), index=cells.index)
    pts = cells[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        inside = MplPath(xy).contains_points(pts)
        out.values[inside & ~out.notna().values] = name
    return out


def process_slide(slide_dir, ann_csv, label):
    """Build a V1-ROI-restricted AnnData for one slide, keeping centroids."""
    cells = pd.read_csv(slide_dir / 'cells.csv.gz').set_index('cell_id')
    polys = load_polygons(ann_csv)
    cells['roi'] = assign(cells, polys)
    cells['group'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None,
        axis=1)
    cells['slide'] = label
    cells['cell_id_orig'] = cells.index
    cells_in = cells.dropna(subset=['group']).copy()

    adata = sc.read_10x_h5(str(slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    adata.obs_names = [f'{label}::{x}' for x in adata.obs_names]
    cells_in.index = [f'{label}::{x}' for x in cells_in.index]
    common = adata.obs_names.intersection(cells_in.index)
    adata = adata[common].copy()
    keep = ['roi', 'group', 'replicate', 'sample_id', 'slide',
            'cell_id_orig', 'x_centroid', 'y_centroid']
    adata.obs = cells_in.loc[common, keep].copy()
    print(f'  [{label}] {adata.n_obs} cells in V1 ROIs')
    return adata


def build_v1(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([a, b], join='outer', merge='same')
    adata.obs_names_make_unique()
    print(f'[build] combined V1: {adata.n_obs} cells x {adata.n_vars} genes')
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slide-a-dir', required=True)
    ap.add_argument('--slide-b-dir', required=True)
    ap.add_argument('--v1-slide-a-ann', required=True)
    ap.add_argument('--v1-slide-b-ann', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--date', default='2026-06-09')
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    D = args.date
    sc.settings.figdir = str(out)
    sc.settings.verbosity = 1

    # ---- build + restrict to V1 ROIs ----
    adata = build_v1(Path(args.slide_a_dir), Path(args.v1_slide_a_ann),
                     Path(args.slide_b_dir), Path(args.v1_slide_b_ann))
    adata.layers['counts'] = adata.X.copy()  # keep raw counts

    # ---- QC (10x workshop) ----
    total = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs['total_counts'] = total
    thres = float(np.quantile(total, MAX_COUNTS_Q))
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
    sc.pp.filter_cells(adata, max_counts=thres)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    print(f'[qc] cells {n0} -> {adata.n_obs} '
          f'(min_counts={MIN_COUNTS}, max_counts={thres:.0f}=q{MAX_COUNTS_Q}); '
          f'genes -> {adata.n_vars} (min_cells={MIN_CELLS})')

    # ---- HVG on raw counts (seurat_v3; skmisc-free reimplementation) ----
    seurat_v3_hvg(adata, n_top_genes=N_HVG, layer='counts')
    print(f'[hvg] {int(adata.var["highly_variable"].sum())} highly-variable '
          f'genes (seurat_v3 via lowess)')

    # ---- normalize / log / scale ----
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()
    sc.pp.scale(adata, zero_center=False, max_value=10)

    # ---- dimensionality reduction + clustering ----
    print('[cluster] pca -> neighbors -> leiden -> umap ...')
    sc.pp.pca(adata, n_comps=N_PCS)
    sc.pp.neighbors(adata, metric='cosine')
    sc.tl.leiden(adata, flavor='igraph', n_iterations=-1, resolution=LEIDEN_RES)
    sc.tl.umap(adata)
    n_clusters = adata.obs['leiden'].nunique()
    print(f'[cluster] {n_clusters} leiden clusters '
          f'(resolution={LEIDEN_RES})')

    # UMAP coloured by cluster (and by slide, to expose any batch effect)
    fig = sc.pl.umap(adata, color=['leiden', 'slide'], show=False,
                     wspace=0.35, return_fig=True)
    fig.savefig(out / f'WholeCell_V1_10x_UMAP_leiden_slide_{D}.pdf',
                bbox_inches='tight')
    plt.close(fig)

    # ---- marker genes per cluster ----
    print('[markers] rank_genes_groups on leiden ...')
    sc.tl.rank_genes_groups(adata, groupby='leiden', layer='lognorm', pts=True)
    fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False,
                                          return_fig=True)
    try:
        fig.savefig(out / f'WholeCell_V1_10x_markers_dotplot_{D}.pdf',
                    bbox_inches='tight')
    except AttributeError:
        plt.savefig(out / f'WholeCell_V1_10x_markers_dotplot_{D}.pdf',
                    bbox_inches='tight')
    plt.close('all')
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(out / f'WholeCell_V1_10x_marker_genes_{D}.csv', index=False)

    # ---- spatial neighbourhood enrichment (per slide) ----
    print('[spatial] spatial_neighbors (delaunay, per slide) -> nhood_enrichment ...')
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
    adata.obs['slide'] = adata.obs['slide'].astype('category')
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True,
                            library_key='slide')
    sq.gr.nhood_enrichment(adata, cluster_key='leiden', seed=0)
    fig, ax = plt.subplots(figsize=(8, 7))
    sq.pl.nhood_enrichment(adata, cluster_key='leiden', cmap='seismic',
                           ax=ax, title='V1 neighbourhood enrichment (leiden)')
    fig.savefig(out / f'WholeCell_V1_10x_nhood_enrichment_{D}.pdf',
                bbox_inches='tight')
    plt.close(fig)

    # ---- exports (write the processed object before the final plot, so the
    #      expensive compute is preserved even if plotting fails) ----
    clustering_res = adata.obs[['cell_id_orig', 'slide', 'leiden']].copy()
    clustering_res.columns = ['cell_id', 'slide', 'group']
    clustering_res.to_csv(out / f'WholeCell_V1_10x_clustering_res_{D}.csv',
                          index=False)
    adata.write(out / f'WholeCell_V1_10x_processed_{D}.h5ad')

    # spatial map of clusters, one panel per slide
    slides = list(adata.obs['slide'].cat.categories)
    fig, axes = plt.subplots(1, len(slides), figsize=(7 * len(slides), 6.5),
                             squeeze=False)
    cats = list(adata.obs['leiden'].cat.categories)
    cmap = plt.get_cmap('tab20', len(cats))
    colmap = {str(c): cmap(i) for i, c in enumerate(cats)}
    leiden_str = adata.obs['leiden'].astype(str)
    for j, sl in enumerate(slides):
        ax = axes[0][j]
        mask = (adata.obs['slide'] == sl).values
        cols = np.array([colmap[c] for c in leiden_str[mask]])
        ax.scatter(adata.obs['x_centroid'].values[mask],
                   adata.obs['y_centroid'].values[mask], s=1.5,
                   c=cols, edgecolor='none', rasterized=True)
        ax.set_title(f'{sl} — leiden clusters', fontsize=11)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])
    handles = [plt.Line2D([0], [0], marker='o', ls='', color=colmap[str(c)],
                          label=str(c), ms=6) for c in cats]
    fig.legend(handles=handles, loc='center right', fontsize=8,
               title='leiden', frameon=False)
    plt.tight_layout(rect=[0, 0, 0.93, 1])
    fig.savefig(out / f'WholeCell_V1_10x_spatial_clusters_{D}.pdf',
                bbox_inches='tight', dpi=200)
    plt.close(fig)

    print(f'[done] outputs in {out}')
    print(f'  cells={adata.n_obs}, genes={adata.n_vars}, '
          f'leiden clusters={n_clusters}')


if __name__ == '__main__':
    main()
