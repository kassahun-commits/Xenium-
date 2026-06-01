#!/usr/bin/env python3
"""
Volcano-plot grid: 6 comparisons × 3 cell types (Excitatory / Inhibitory /
Astrocyte) = 18 volcanos.

Comparisons (per cell type):
  1. EtOH_veh         vs H2O_veh        (acute alcohol)
  2. ChronicEtOH      vs H2O_veh        (chronic alcohol)
  3. H2O_MCT1i        vs H2O_veh        (drug-only)
  4. EtOH_MCT1i       vs H2O_veh        (alcohol+drug vs control)
  5. EtOH_MCT1i       vs EtOH_veh       (MCT1i rescue of acute alcohol)
  6. MAT2A_OE         vs MAT2A_CM       (MAT2A OE vs catalytic-mutant control)

Layout: 1 comparison per page; 3 cell-type panels per page (E | I | A).
Stats: cell-level Wilcoxon (scanpy.tl.rank_genes_groups). Same
pseudoreplication caveat as before.

Outputs:
  - SlidesAB_volcanos_AllContrasts_<date>.pdf  (multi-page)
  - SlidesAB_volcanos_DE_combined_<date>.csv   (long-form DE table)
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
    # APP added to V1 Slide A
    'APP 1':          ('APP',        1),
    'APP 2':          ('APP',        2),
    'APP 3':          ('APP',        3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']
INHIB_MARKERS = ['Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln', 'Lhx6']

SUBTYPES = ['Excitatory', 'Inhibitory', 'Astrocyte']
SUBTYPE_COLORS = {'Excitatory': '#1f77b4',
                  'Inhibitory': '#9467bd',
                  'Astrocyte':  '#d62728'}

# (test_group, reference_group, short_label)
COMPARISONS = [
    ('EtOH_veh',    'H2O_veh',  'Acute alcohol — EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic alcohol — ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'Drug-only — H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'Alcohol + drug — EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'MCT1i rescue — EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2A overexpression — MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP — APP vs H2O_veh'),
]

# Significance thresholds (defaults; overridable via CLI)
LFC_THRESH = 1.0
PADJ_THRESH = 1e-3


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


def run_de(adata_sub, test, ref):
    """Run cell-level Wilcoxon; return DataFrame and meta."""
    d = adata_sub[adata_sub.obs['group'].astype(str).isin([test, ref])].copy()
    if d.shape[0] < 30:
        return None, {'n_test_cells': 0, 'n_ref_cells': 0, 'n_test_mice': 0, 'n_ref_mice': 0}
    d.obs['de_group'] = d.obs['group'].astype(str)
    n_t = int((d.obs['de_group'] == test).sum())
    n_r = int((d.obs['de_group'] == ref).sum())
    s_t = sorted(d.obs.loc[d.obs['de_group'] == test, 'sample_id'].unique().tolist())
    s_r = sorted(d.obs.loc[d.obs['de_group'] == ref,  'sample_id'].unique().tolist())
    if min(n_t, n_r) < 10:
        return None, {'n_test_cells': n_t, 'n_ref_cells': n_r,
                      'n_test_mice': len(s_t), 'n_ref_mice': len(s_r)}
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
    return df, {'n_test_cells': n_t, 'n_ref_cells': n_r,
                'n_test_mice': len(s_t), 'n_ref_mice': len(s_r)}


def plot_volcano(ax, de_df, title, meta, subtitle_color):
    """Render one volcano on the given axis."""
    if de_df is None or len(de_df) == 0:
        ax.text(0.5, 0.5, '(insufficient cells)', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=10, color=subtitle_color)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ('top','right'): ax.spines[s].set_visible(False)
        return

    x = de_df['logfc'].values
    y = -np.log10(np.clip(de_df['padj'].values, 1e-300, 1))
    sig = (de_df['logfc'].abs() >= LFC_THRESH) & (de_df['padj'] < PADJ_THRESH)

    # Background (non-sig)
    ax.scatter(x[~sig], y[~sig], s=4, c='lightgray', alpha=0.45, edgecolor='none')
    # Up in test (red), up in ref (blue)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=12, c='#d62728',
               alpha=0.8, edgecolor='none')
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=12, c='#1f77b4',
               alpha=0.8, edgecolor='none')
    # Labels for top 12 each direction
    top_up = de_df[sig & (de_df['logfc'] > 0)].nlargest(12, 'logfc')
    top_dn = de_df[sig & (de_df['logfc'] < 0)].nsmallest(12, 'logfc')
    for _, row in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(row['gene'], (row['logfc'], -np.log10(max(row['padj'], 1e-300))),
                    fontsize=6, alpha=0.85)
    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(LFC_THRESH,  ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(-LFC_THRESH, ls='--', color='black', lw=0.4, alpha=0.4)
    ax.set_xlabel('log2 fold change', fontsize=8)
    ax.set_ylabel('-log10 adj p-value', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    n_up = int((sig & (x > 0)).sum()); n_dn = int((sig & (x < 0)).sum())
    n_info = (f"test cells={meta['n_test_cells']:,} / mice={meta['n_test_mice']}\n"
              f"ref cells={meta['n_ref_cells']:,} / mice={meta['n_ref_mice']}\n"
              f"up={n_up}, down={n_dn}")
    ax.text(0.02, 0.98, n_info, transform=ax.transAxes, fontsize=6,
            va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='gray', alpha=0.85))
    ax.set_title(title, fontsize=10, color=subtitle_color, fontweight='bold')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--lfc-thresh', type=float, default=1.0,
                    help='|log2 fold change| threshold for significance (default 1.0)')
    ap.add_argument('--padj-thresh', type=float, default=1e-3,
                    help='Adjusted p-value threshold (default 1e-3)')
    ap.add_argument('--suffix', default='',
                    help='Optional suffix added to output filenames')
    args = ap.parse_args()
    global LFC_THRESH, PADJ_THRESH
    LFC_THRESH = args.lfc_thresh
    PADJ_THRESH = args.padj_thresh
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_volcanos_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Volcano grid: 6 comparisons x 3 cell types ===')

    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Cell-typing
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

    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory',
                      random_state=0, n_bins=25)
    inhib_present = [g for g in INHIB_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=inhib_present, score_name='score_Inhibitory',
                      random_state=0, n_bins=25)
    is_excit = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Excitatory'] > 0))
    is_inhib = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Inhibitory'] > 0))
    is_astro = (adata.obs['celltype'].astype(str) == 'Astrocyte')

    subtype_adatas = {
        'Excitatory': adata[is_excit].copy(),
        'Inhibitory': adata[is_inhib].copy(),
        'Astrocyte':  adata[is_astro].copy(),
    }
    for st, a in subtype_adatas.items():
        logging.info('%s subset: %d cells', st, a.shape[0])

    # Run all DEs
    de_long_rows = []
    de_cache = {}  # (subtype, test, ref) -> (df, meta)
    for subtype, a in subtype_adatas.items():
        for test, ref, label in COMPARISONS:
            df, meta = run_de(a, test, ref)
            de_cache[(subtype, test, ref)] = (df, meta)
            if df is not None:
                df_out = df.copy()
                df_out['subtype'] = subtype
                df_out['test'] = test; df_out['reference'] = ref
                df_out['n_test_cells']  = meta['n_test_cells']
                df_out['n_ref_cells']   = meta['n_ref_cells']
                df_out['n_test_mice']   = meta['n_test_mice']
                df_out['n_ref_mice']    = meta['n_ref_mice']
                de_long_rows.append(df_out)
            logging.info('DE %s: %s vs %s -> %s',
                         subtype, test, ref,
                         'OK' if df is not None else f"skipped (cells {meta['n_test_cells']}/{meta['n_ref_cells']})")
    if de_long_rows:
        de_long = pd.concat(de_long_rows, ignore_index=True)
        de_csv = args.out_dir / f'{args.label}_volcanos_DE_combined_{args.date}{args.suffix}.csv'
        de_long.to_csv(de_csv, index=False)
        logging.info('Wrote %s', de_csv)

    # ---- Plot ----
    pdf_out = args.out_dir / f'{args.label}_volcanos_AllContrasts_{args.date}{args.suffix}.pdf'
    with PdfPages(pdf_out) as pdf:
        for test, ref, label in COMPARISONS:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
            for ax, subtype in zip(axes, SUBTYPES):
                df, meta = de_cache.get((subtype, test, ref), (None, {}))
                title = f'{subtype}'
                plot_volcano(ax, df, title, meta, SUBTYPE_COLORS[subtype])
            fig.suptitle(label + '\n(red = up in test, blue = up in reference; '
                                  f'sig: |log2FC|>={LFC_THRESH}, padj<{PADJ_THRESH:g})',
                          fontsize=11, y=1.02)
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    logging.info('Wrote %s', pdf_out)


if __name__ == '__main__':
    main()
