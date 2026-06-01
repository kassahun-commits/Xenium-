#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled (V1) — ENDOTHELIAL-cell DE heat map.

Endothelial cells are pulled out of the "Other" bucket by adding a pan-
endothelial marker module to the argmax cell-typing step, then a cell-level
Wilcoxon DE is run for five contrasts and the top differentially expressed
genes are shown as a log2 fold-change heat map (rows = genes, columns =
contrasts).

Contrasts (in column order, as requested):
  1. MCT1i        = H2O_MCT1i   vs H2O_veh   (drug only)
  2. EtOH         = EtOH_veh    vs H2O_veh   (acute alcohol)
  3. EtOH+MCT1i   = EtOH_MCT1i  vs H2O_veh   (alcohol + drug vs control)
  4. EtOH+MCT1i / EtOH = EtOH_MCT1i vs EtOH_veh (MCT1i on top of alcohol)
  5. Chronic EtOH = ChronicEtOH vs H2O_veh   (chronic alcohol)

Row (gene) selection: union, over the five contrasts, of the top TOP_N genes
ranked by adjusted p-value among genes passing |log2FC| >= LFC_THRESH and
padj < SEL_PADJ.  Rows ordered by the acute-EtOH (EtOH_veh vs H2O_veh) log2FC,
strongest up on top.  Cell colour = log2FC (RdBu_r, clipped to +/-CLIP);
significance asterisks (* padj<0.05, ** <0.01, *** <0.001) per cell.

Stats: cell-level Wilcoxon (scanpy.tl.rank_genes_groups, BH-corrected).  Same
pseudoreplication caveat as the rest of the project — treat direction and
log2FC magnitude as primary; adjusted p-values are anticonservative.  NOTE:
endothelial cells are a minority population, so per-group counts are logged and
any contrast with too few cells is left blank (grey).

V1 only.  All paths via CLI.  Vector PDF, editable text (pdf.fonttype=42); one
source-data CSV (log2FC matrix) + one long-form DE table written alongside.
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
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# --------------------------------------------------------------------------- #
#  parameters                                                                 #
# --------------------------------------------------------------------------- #
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
    'APP 1':          ('APP',        1),
    'APP 2':          ('APP',        2),
    'APP 3':          ('APP',        3),
}

# cell-typing marker modules; argmax over these assigns a coarse type.
# Endothelial added so vascular cells stop falling into "Other".
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
    'Endothelial':     ['Pecam1', 'Cdh5', 'Flt1', 'Kdr', 'Tek', 'Vwf',
                        'Sox17', 'Cd34', 'Eng', 'Mfsd2a'],
}

# (test, reference, column label) — column order as requested by user
CONTRASTS = [
    ('H2O_MCT1i',   'H2O_veh',  'MCT1i\nvs H2O'),
    ('EtOH_veh',    'H2O_veh',  'EtOH\nvs H2O'),
    ('EtOH_MCT1i',  'H2O_veh',  'EtOH+MCT1i\nvs H2O'),
    ('EtOH_MCT1i',  'EtOH_veh', 'EtOH+MCT1i\nvs EtOH'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic EtOH\nvs H2O'),
]
ORDER_TEST, ORDER_REF = 'EtOH_veh', 'H2O_veh'   # rows ordered by this contrast's log2FC

LFC_THRESH = 0.5      # |log2FC| for a gene to count as DE
SEL_PADJ   = 0.05     # adj-p for row selection (looser than 1e-3: endo cells are rare)
TOP_N      = 10       # top genes per contrast (by padj) -> unioned across contrasts
CLIP       = 3.0      # colour-scale clip (log2FC)
MIN_CELLS  = 10       # per-group minimum for a contrast to be run


# --------------------------------------------------------------------------- #
#  loaders / typing (shared with the rest of the project)                     #
# --------------------------------------------------------------------------- #
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
    cells['group'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)
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


def type_cells(adata):
    """argmax marker-module typing; returns coarse 'celltype' on adata.obs."""
    for ct, gs in MARKERS.items():
        present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=present, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
        missing = [g for g in gs if g not in adata.var_names]
        if missing:
            logging.info('  %s module: %d/%d markers present (missing: %s)',
                         ct, len(present), len(gs), ','.join(missing))
    score_cols = [f'score_{ct}' for ct in MARKERS]
    scores = adata.obs[score_cols]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_val = scores.max(axis=1)
    best_ct[best_val <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(
        best_ct, categories=list(MARKERS.keys()) + ['Unclassified'])
    return adata


def run_de(adata_sub, test, ref):
    """Cell-level Wilcoxon test vs ref. Returns (DataFrame|None, meta)."""
    d = adata_sub[adata_sub.obs['group'].astype(str).isin([test, ref])].copy()
    n_t = int((d.obs['group'].astype(str) == test).sum())
    n_r = int((d.obs['group'].astype(str) == ref).sum())
    s_t = sorted(d.obs.loc[d.obs['group'].astype(str) == test, 'sample_id'].unique().tolist())
    s_r = sorted(d.obs.loc[d.obs['group'].astype(str) == ref,  'sample_id'].unique().tolist())
    meta = {'n_test_cells': n_t, 'n_ref_cells': n_r,
            'n_test_mice': len(s_t), 'n_ref_mice': len(s_r)}
    if min(n_t, n_r) < MIN_CELLS or (n_t + n_r) < 30:
        return None, meta
    d.obs['de_group'] = d.obs['group'].astype(str)
    sc.tl.rank_genes_groups(d, 'de_group', method='wilcoxon',
                            reference=ref, n_genes=d.shape[1])
    r = d.uns['rank_genes_groups']
    df = pd.DataFrame({
        'gene':  [str(g) for g in r['names'][test]],
        'logfc': r['logfoldchanges'][test],
        'pval':  r['pvals'][test],
        'padj':  r['pvals_adj'][test],
        'score': r['scores'][test],
    })
    return df, meta


# --------------------------------------------------------------------------- #
#  gene selection + plotting                                                  #
# --------------------------------------------------------------------------- #
def select_genes(de_by_label):
    """Union of top-N (by padj) DE genes across contrasts; ordered by acute-EtOH log2FC."""
    selected = set()
    for _, _, lab in CONTRASTS:
        df = de_by_label.get(lab)
        if df is None:
            continue
        sig = df[(df['logfc'].abs() >= LFC_THRESH) & (df['padj'] < SEL_PADJ)]
        top = sig.sort_values('padj').head(TOP_N)
        selected.update(top['gene'].tolist())
    if not selected:
        return [], None, None
    # ordering key = acute-EtOH log2FC (fallback: mean over available contrasts)
    order_label = next(l for t, r, l in CONTRASTS if (t, r) == (ORDER_TEST, ORDER_REF))
    order_df = de_by_label.get(order_label)
    genes = list(selected)
    if order_df is not None:
        key = order_df.set_index('gene')['logfc'].reindex(genes)
    else:
        key = pd.Series(0.0, index=genes)
    key = key.fillna(0.0).sort_values(ascending=False)
    ordered = list(key.index)
    # build log2FC + padj matrices (genes x contrasts)
    lfc_cols, padj_cols = {}, {}
    for _, _, lab in CONTRASTS:
        df = de_by_label.get(lab)
        if df is None:
            lfc_cols[lab] = pd.Series(np.nan, index=ordered)
            padj_cols[lab] = pd.Series(np.nan, index=ordered)
        else:
            s = df.set_index('gene')
            lfc_cols[lab] = s['logfc'].reindex(ordered)
            padj_cols[lab] = s['padj'].reindex(ordered)
    lfc_mat = pd.DataFrame(lfc_cols, index=ordered)
    padj_mat = pd.DataFrame(padj_cols, index=ordered)
    return ordered, lfc_mat, padj_mat


def stars(p):
    if pd.isna(p):
        return ''
    if p < 1e-3:
        return '***'
    if p < 1e-2:
        return '**'
    if p < 5e-2:
        return '*'
    return ''


def draw_page(pdf, lfc_mat, padj_mat, counts, today):
    col_labels = [lab for _, _, lab in CONTRASTS]
    genes = list(lfc_mat.index)
    ng = len(genes)
    fig_h = max(5.0, 0.20 * ng + 2.4)
    fig, ax = plt.subplots(figsize=(6.4, fig_h))
    vals = lfc_mat.values.astype(float)
    masked = np.ma.masked_invalid(vals)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('lightgray')
    im = ax.imshow(masked, aspect='auto', cmap=cmap, vmin=-CLIP, vmax=CLIP)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticks(range(ng))
    ax.set_yticklabels(genes, fontsize=max(4, min(8, 460 / max(ng, 1))), style='italic')
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # significance asterisks
    for i in range(ng):
        for j in range(len(col_labels)):
            mark = stars(padj_mat.values[i, j])
            if mark:
                v = vals[i, j]
                txt_c = 'white' if (not np.isnan(v) and abs(v) > CLIP * 0.55) else 'black'
                ax.text(j, i, mark, ha='center', va='center',
                        fontsize=5.5, color=txt_c)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03, extend='both')
    cbar.set_label('log2 fold change', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    count_str = '   '.join(f'{g}={counts.get(g, 0):,}' for g in
                           ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i', 'ChronicEtOH'])
    fig.suptitle('V1 — Endothelial cells — top DE genes across contrasts\n'
                 f'log2FC vs reference; rows = union of top {TOP_N}/contrast '
                 f'(|log2FC|>={LFC_THRESH:g}, padj<{SEL_PADJ:g}); colour clipped +/-{CLIP:g}',
                 fontsize=9.5, y=0.995)
    ax.text(0.0, -0.085 - 0.012 * (ng > 30),
            f'endothelial cells/group:  {count_str}\n'
            f'* padj<0.05  ** <0.01  *** <0.001   |   cell-level Wilcoxon (BH); '
            f'grey = too few cells.  {today} · MEWS Lab',
            transform=ax.transAxes, fontsize=6, va='top', color='#444444')
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  main                                                                       #
# --------------------------------------------------------------------------- #
def main():
    global CLIP, TOP_N
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--clip', type=float, default=CLIP)
    ap.add_argument('--top-n', type=int, default=TOP_N)
    ap.add_argument('--date', default=date.today().isoformat())
    args = ap.parse_args()

    CLIP = args.clip
    TOP_N = args.top_n
    args.outdir.mkdir(parents=True, exist_ok=True)

    log_path = args.outdir / f'EndothelialDE_Heatmap_V1_log_{args.date}.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Endothelial DE heat map (V1, Slides A+B pooled) ===')

    ada_a = process_slide(args.slide_a_dir, args.v1_slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.v1_slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    logging.info('Pooled V1: %d cells x %d genes', adata.shape[0], adata.shape[1])

    type_cells(adata)
    comp = adata.obs['celltype'].value_counts().to_dict()
    logging.info('Cell-type composition: %s', comp)

    endo = adata[adata.obs['celltype'].astype(str) == 'Endothelial'].copy()
    counts = endo.obs['group'].astype(str).value_counts().to_dict()
    logging.info('Endothelial cells: %d total', endo.shape[0])
    logging.info('Endothelial per group: %s', counts)

    # run DE per contrast
    de_by_label, de_long_rows = {}, []
    for test, ref, lab in CONTRASTS:
        df, meta = run_de(endo, test, ref)
        de_by_label[lab] = df
        if df is None:
            logging.warning('Contrast %s vs %s SKIPPED (test=%d ref=%d cells)',
                            test, ref, meta['n_test_cells'], meta['n_ref_cells'])
            continue
        logging.info('Contrast %s vs %s: test=%d ref=%d cells', test, ref,
                     meta['n_test_cells'], meta['n_ref_cells'])
        out = df.copy()
        out['cell_type'] = 'Endothelial'
        out['test'] = test
        out['reference'] = ref
        out['contrast'] = lab.replace('\n', ' ')
        out['n_test_cells'] = meta['n_test_cells']
        out['n_ref_cells'] = meta['n_ref_cells']
        out['n_test_mice'] = meta['n_test_mice']
        out['n_ref_mice'] = meta['n_ref_mice']
        de_long_rows.append(out)

    if de_long_rows:
        de_long = pd.concat(de_long_rows, ignore_index=True)
        de_csv = args.outdir / f'SlidesAB_EndothelialDE_5contrasts_V1_{args.date}.csv'
        de_long.to_csv(de_csv, index=False)
        logging.info('Wrote %s (%d rows)', de_csv, len(de_long))

    ordered, lfc_mat, padj_mat = select_genes(de_by_label)
    if not ordered:
        logging.error('No genes passed selection thresholds — nothing to plot.')
        return
    logging.info('Selected %d genes for heat map', len(ordered))

    # source-data CSV (log2FC matrix + padj)
    src = lfc_mat.copy()
    src.columns = [c.replace('\n', ' ') + ' (log2FC)' for c in src.columns]
    for (_, _, lab) in CONTRASTS:
        src[lab.replace('\n', ' ') + ' (padj)'] = padj_mat[lab].values
    src.insert(0, 'gene', src.index)
    src_csv = args.outdir / f'SlidesAB_EndothelialDE_Heatmap_5contrasts_V1_{args.date}.csv'
    src.to_csv(src_csv, index=False)
    logging.info('Wrote %s', src_csv)

    pdf_out = args.outdir / f'Xenium_EndothelialDE_Heatmap_5contrasts_V1_{args.date}.pdf'
    with PdfPages(pdf_out) as pdf:
        draw_page(pdf, lfc_mat, padj_mat, counts, args.date)
    logging.info('Wrote %s', pdf_out)
    logging.info('Done.')


if __name__ == '__main__':
    main()
