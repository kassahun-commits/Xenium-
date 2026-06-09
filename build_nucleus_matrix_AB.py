#!/usr/bin/env python3
"""
Xenium_May2026 — Nucleus-only count matrix for Slides A + B.

Purpose
-------
Produce a version of the data restricted to *nuclear* transcripts, so it can be
compared against the standard whole-cell `cell_feature_matrix.h5` (which 10x
builds from ALL transcripts assigned to a cell — nucleus + cytoplasm).

Why this lives at the transcript level
--------------------------------------
`overlaps_nucleus` is a per-transcript flag in `transcripts.parquet`, not a
cell-level attribute. So we rebuild the cell x gene count matrix from scratch,
counting only transcripts that fall inside the nucleus, then feed the result
through the SAME downstream steps as the existing whole-cell pipeline
(slideB_alcohol_prelim.py / slidesAB_combined_curated_dotplot.py) so the two are
directly comparable.

Pipeline
--------
  1. Lock the gene set to the whole-cell matrix var_names (comparability).
  2. Stream transcripts.parquet by row-group; keep a transcript only if:
       overlaps_nucleus == 1   (in nucleus)
       qv >= --qv-min          (default 20; data-driven, see notes below)
       is_gene == True         (drop negative-control / unassigned codewords)
       cell_id != 'UNASSIGNED' (assigned to a real cell)
     Accumulate per-(cell, gene) counts -> sparse nucleus matrix.
  3. ROI assignment: point-in-polygon on cells.csv.gz centroids using the
     Xenium Explorer annotation polygons; map ROI label -> (group, replicate)
     via LABEL_MAP. Unmapped ROIs (e.g. APP 1/2/3 on Slide A) are excluded,
     matching the existing pipeline.
  4. Concatenate Slide A + Slide B, QC filter, normalize to 1e4, log1p.
  5. Outputs: per-slide + combined .h5ad (raw counts kept in layers['counts']),
     per-cell ROI tables (.csv), a QC verification figure (vector PDF, editable
     text) with source .csv, and a full log.

QC parameter rationale (--qv-min default 20)
---------------------------------------------
The Slide A QV distribution (gene transcripts) is bimodal: a low-confidence
"noise" mode peaking < 10, a high-confidence "signal" mode piling up at 35-40,
and a trough at QV 18-20. qv >= 20 sits at that natural valley AND matches the
cutoff 10x used to build the whole-cell matrix, keeping the nucleus matrix
comparable. min_counts=10 / min_cells=5 match the established pipeline.

Paths are passed via CLI (MEWS Lab rule 4: no hardcoded paths).
Resumable: per-slide nucleus counts are cached to .npz; re-running skips the
expensive parquet pass if the cache exists (MEWS Lab rule 7).
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
import pyarrow.parquet as pq
import scanpy as sc
import scipy.sparse as sp
from matplotlib.path import Path as MplPath

# Editable text in vector outputs (MEWS Lab rule 4.3)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

UNASSIGNED = 'UNASSIGNED'
TRANSCRIPT_COLS = ['cell_id', 'feature_name', 'qv', 'is_gene', 'overlaps_nucleus']

# Translate Xenium Explorer labels -> (group, replicate). Both slides.
# (Copied verbatim from slidesAB_combined_curated_dotplot.py for consistency.)
LABEL_MAP = {
    # --- Slide B labels ---
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
    # --- Slide A labels ---
    'EtOH +MCT1i 2':  ('EtOH_MCT1i', 2),
    'EtOH +MCT1i 3':  ('EtOH_MCT1i', 3),
    'H20+veh 2':      ('H2O_veh',    2),
    'H20 +MCT1i 2':   ('H2O_MCT1i',  2),
    'Chronic EtOH 2': ('ChronicEtOH', 2),
    'chronic EtOH 3': ('ChronicEtOH', 3),
}

GROUP_ORDER = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
               'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']


# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--slide-a-dir', required=True, type=Path,
                   help='Xenium output folder for Slide A (0088372)')
    p.add_argument('--slide-a-ann', required=True, type=Path,
                   help='Slide A polygon annotations CSV')
    p.add_argument('--slide-b-dir', required=True, type=Path,
                   help='Xenium output folder for Slide B (0088233)')
    p.add_argument('--slide-b-ann', required=True, type=Path,
                   help='Slide B polygon annotations CSV')
    p.add_argument('--out-dir', required=True, type=Path)
    p.add_argument('--label', default='SlidesAB_nucleus')
    p.add_argument('--qv-min', type=float, default=20.0,
                   help='Min Q-Score per transcript (default 20; data-driven)')
    p.add_argument('--min-counts', type=int, default=10,
                   help='Min nuclear transcripts per cell (QC)')
    p.add_argument('--min-cells', type=int, default=5,
                   help='Min cells expressing a gene (QC)')
    p.add_argument('--norm-target', type=float, default=1e4)
    p.add_argument('--date', default=date.today().isoformat())
    p.add_argument('--force', action='store_true',
                   help='Ignore cached per-slide counts and rebuild')
    return p.parse_args()


def setup_logging(out_dir: Path, label: str) -> Path:
    log_path = out_dir / f'{label}_log_{date.today().isoformat()}.txt'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
        handlers=[logging.FileHandler(log_path, mode='w'),
                  logging.StreamHandler(sys.stdout)])
    return log_path


def get_gene_index(slide_dir: Path) -> tuple[list[str], dict[str, int]]:
    """Gene set = var_names of the whole-cell matrix (ensures comparability)."""
    adata = sc.read_10x_h5(str(slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    genes = adata.var_names.tolist()
    return genes, {g: i for i, g in enumerate(genes)}


def build_nucleus_counts(slide_dir: Path, gene_to_idx: dict[str, int],
                         qv_min: float, cache_npz: Path,
                         force: bool = False):
    """Stream transcripts.parquet -> nucleus-only cell x gene CSR matrix.

    Returns (csr_matrix [n_cells x n_genes], list_of_cell_ids, stats dict).
    Resumable: result cached to `cache_npz` (+ sibling .cells.npy).
    """
    cells_npy = cache_npz.with_suffix('.cells.npy')
    if cache_npz.exists() and cells_npy.exists() and not force:
        logging.info('  [cache] loading nucleus counts from %s', cache_npz)
        mat = sp.load_npz(cache_npz)
        cell_ids = np.load(cells_npy, allow_pickle=True).tolist()
        return mat, cell_ids, {'cached': True}

    n_genes = len(gene_to_idx)
    path = slide_dir / 'transcripts.parquet'
    pf = pq.ParquetFile(path)
    n_rg = pf.num_row_groups

    cell_to_idx: dict[str, int] = {}
    key_chunks: list[np.ndarray] = []
    cnt_chunks: list[np.ndarray] = []
    n_seen = n_kept = 0

    def consolidate():
        """Collapse accumulated (key,count) chunks to keep memory bounded."""
        if not key_chunks:
            return
        allk = np.concatenate(key_chunks)
        allc = np.concatenate(cnt_chunks)
        uk, inv = np.unique(allk, return_inverse=True)
        agg = np.zeros(len(uk), dtype=np.int64)
        np.add.at(agg, inv, allc)
        key_chunks.clear(); cnt_chunks.clear()
        key_chunks.append(uk); cnt_chunks.append(agg)

    for rg in range(n_rg):
        df = pf.read_row_group(rg, columns=TRANSCRIPT_COLS).to_pandas()
        n_seen += len(df)
        m = (df['overlaps_nucleus'].to_numpy() == 1) & \
            df['is_gene'].to_numpy() & \
            (df['qv'].to_numpy() >= qv_min) & \
            (df['cell_id'].to_numpy() != UNASSIGNED)
        df = df.loc[m, ['cell_id', 'feature_name']]
        if df.empty:
            continue
        gi = df['feature_name'].map(gene_to_idx).to_numpy()
        ok = ~pd.isna(gi)              # drop genes not in the locked panel
        if not ok.all():
            df = df[ok]; gi = gi[ok]
        gi = gi.astype(np.int64)
        n_kept += len(df)

        cell_arr = df['cell_id'].to_numpy()
        uniq, inv = np.unique(cell_arr, return_inverse=True)
        glob = np.empty(len(uniq), dtype=np.int64)
        for i, c in enumerate(uniq):
            j = cell_to_idx.get(c)
            if j is None:
                j = len(cell_to_idx); cell_to_idx[c] = j
            glob[i] = j
        cell_idx = glob[inv]

        key = cell_idx * n_genes + gi
        uk, cnt = np.unique(key, return_counts=True)
        key_chunks.append(uk); cnt_chunks.append(cnt.astype(np.int64))

        if rg % 10 == 9:
            consolidate()
            logging.info('  row-group %d/%d | seen=%s kept=%s cells=%s',
                         rg + 1, n_rg, f'{n_seen:,}', f'{n_kept:,}',
                         f'{len(cell_to_idx):,}')

    consolidate()
    n_cells = len(cell_to_idx)
    if not key_chunks or n_cells == 0:
        raise RuntimeError('No nuclear transcripts passed filters.')
    keys = key_chunks[0]; data = cnt_chunks[0].astype(np.float32)
    rows = keys // n_genes
    cols = keys % n_genes
    mat = sp.coo_matrix((data, (rows, cols)),
                        shape=(n_cells, n_genes)).tocsr()
    cell_ids = [None] * n_cells
    for c, j in cell_to_idx.items():
        cell_ids[j] = c

    cache_npz.parent.mkdir(parents=True, exist_ok=True)
    sp.save_npz(cache_npz, mat)
    np.save(cells_npy, np.array(cell_ids, dtype=object))
    stats = {'cached': False, 'n_seen': n_seen, 'n_kept': n_kept,
             'n_cells': n_cells}
    logging.info('  built nucleus matrix: %s cells x %s genes | '
                 'kept %s / %s transcripts (%.1f%%)',
                 f'{n_cells:,}', f'{n_genes:,}', f'{n_kept:,}',
                 f'{n_seen:,}', 100.0 * n_kept / max(n_seen, 1))
    return mat, cell_ids, stats


def load_polygons(csv_path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(csv_path, comment='#')
    if not {'Selection', 'X', 'Y'}.issubset(df.columns):
        raise ValueError(f'Annotations missing cols; got {df.columns.tolist()}')
    return {name: sub[['X', 'Y']].values.astype(float)
            for name, sub in df.groupby('Selection', sort=False)}


def assign_cells_to_rois(cells_df: pd.DataFrame,
                         polygons: dict[str, np.ndarray]) -> pd.Series:
    """Point-in-polygon on centroids. First assignment wins on overlap."""
    out = pd.Series(np.full(len(cells_df), np.nan, dtype=object),
                    index=cells_df.index, name='roi')
    pts = cells_df[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        inside = MplPath(xy).contains_points(pts)
        clash = inside & out.notna().values
        if clash.any():
            logging.warning('  %s cells already assigned also inside "%s" '
                            '- keeping first', int(clash.sum()), name)
            inside = inside & ~out.notna().values
        out.values[inside] = name
    return out


def build_slide_adata(slide_dir: Path, ann_path: Path, slide_label: str,
                      genes: list[str], gene_to_idx: dict[str, int],
                      qv_min: float, out_dir: Path, force: bool) -> ad.AnnData:
    logging.info('=== %s ===', slide_label)
    cache = out_dir / f'_cache_{slide_label}_nucleus_counts_qv{int(qv_min)}.npz'
    mat, cell_ids, _ = build_nucleus_counts(
        slide_dir, gene_to_idx, qv_min, cache, force=force)

    adata = ad.AnnData(X=mat,
                       obs=pd.DataFrame(index=pd.Index(cell_ids, name='cell_id')),
                       var=pd.DataFrame(index=pd.Index(genes, name='gene')))

    # ROI assignment from centroids
    cells = pd.read_csv(slide_dir / 'cells.csv.gz').set_index('cell_id')
    polygons = load_polygons(ann_path)
    logging.info('  polygons: %d | cells with centroids: %s',
                 len(polygons), f'{len(cells):,}')
    cells['roi'] = assign_cells_to_rois(cells, polygons)
    cells['group'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['slide'] = slide_label
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None,
        axis=1)

    unmapped = set(cells['roi'].dropna().unique()) - set(LABEL_MAP.keys())
    if unmapped:
        logging.info('  unmapped ROIs excluded (expected, e.g. APP): %s',
                     sorted(unmapped))

    cells_in = cells.dropna(subset=['group']).copy()
    common = adata.obs_names.intersection(cells_in.index)
    logging.info('  cells with nuclear counts AND a known group: %s',
                 f'{len(common):,}')
    adata = adata[common].copy()
    adata.obs = cells_in.loc[common, ['roi', 'group', 'replicate',
                                      'sample_id', 'slide',
                                      'x_centroid', 'y_centroid']].copy()
    # Prefix to avoid cell_id collisions across slides on concat
    adata.obs_names = [f'{slide_label}::{x}' for x in adata.obs_names]
    return adata


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(out_dir, args.label)
    L, D = args.label, args.date

    logging.info('=== Nucleus-only matrix build: Slides A + B ===')
    logging.info('qv_min=%s | min_counts=%d | min_cells=%d | norm=%g',
                 args.qv_min, args.min_counts, args.min_cells, args.norm_target)

    genes, gene_to_idx = get_gene_index(args.slide_a_dir)
    logging.info('gene panel locked to whole-cell matrix: %d genes', len(genes))

    a = build_slide_adata(args.slide_a_dir, args.slide_a_ann, 'SlideA',
                          genes, gene_to_idx, args.qv_min, out_dir, args.force)
    b = build_slide_adata(args.slide_b_dir, args.slide_b_ann, 'SlideB',
                          genes, gene_to_idx, args.qv_min, out_dir, args.force)

    # Save per-slide (pre-QC) objects
    a.write_h5ad(out_dir / f'{L}_SlideA_preQC_{D}.h5ad')
    b.write_h5ad(out_dir / f'{L}_SlideB_preQC_{D}.h5ad')

    adata = ad.concat([a, b], join='outer', label=None, index_unique=None)
    adata.obs['group'] = pd.Categorical(adata.obs['group'].astype(str),
                                        categories=GROUP_ORDER, ordered=True)
    logging.info('Combined (pre-QC): %s cells x %s genes',
                 f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

    # ---- QC ----
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    n0, g0 = adata.shape
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    logging.info('After QC: %s cells (-%s) x %s genes (-%s)',
                 f'{adata.shape[0]:,}', f'{n0 - adata.shape[0]:,}',
                 f'{adata.shape[1]:,}', f'{g0 - adata.shape[1]:,}')

    # ---- Normalize + log ----
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=args.norm_target)
    sc.pp.log1p(adata)

    # ---- Outputs ----
    h5 = out_dir / f'{L}_combined_{D}.h5ad'
    adata.write_h5ad(h5)
    logging.info('Wrote %s', h5)

    cell_csv = out_dir / f'{L}_combined_cells_with_roi_{D}.csv'
    obs_out = adata.obs.copy()
    obs_out['total_counts_nucleus'] = adata.layers['counts'].sum(axis=1).A1
    obs_out.to_csv(cell_csv)
    logging.info('Wrote %s', cell_csv)

    # Per-sample cell counts (verification)
    samp = (adata.obs.groupby(['slide', 'group', 'replicate'], observed=True)
            .size().rename('n_cells').reset_index().sort_values(['group', 'slide']))
    samp_csv = out_dir / f'{L}_cells_per_sample_{D}.csv'
    samp.to_csv(samp_csv, index=False)
    logging.info('Wrote %s\n%s', samp_csv, samp.to_string(index=False))

    # ---- QC verification figure (vector PDF, editable text) ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    tc = np.asarray(adata.obs['total_counts'])
    axes[0].hist(tc, bins=80, color='#444')
    axes[0].set_xlabel('Nuclear transcripts per cell (post-QC)')
    axes[0].set_ylabel('cells'); axes[0].set_yscale('log')
    axes[0].set_title(f'{L}: counts per cell')
    s = samp.copy()
    s['lab'] = [f'{g}_{int(r)} ({sl})' for g, r, sl in
                zip(s['group'].astype(str), s['replicate'], s['slide'].astype(str))]
    axes[1].bar(range(len(s)), s['n_cells'], color='#3b7', edgecolor='black',
                linewidth=0.4)
    axes[1].set_xticks(range(len(s)))
    axes[1].set_xticklabels(s['lab'], rotation=45, ha='right', fontsize=7)
    axes[1].set_ylabel('cells'); axes[1].set_title(f'{L}: cells per sample')
    plt.tight_layout()
    qc_pdf = out_dir / f'{L}_QCverification_{D}.pdf'
    plt.savefig(qc_pdf, bbox_inches='tight'); plt.close()
    logging.info('Wrote %s', qc_pdf)

    logging.info('Done.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
