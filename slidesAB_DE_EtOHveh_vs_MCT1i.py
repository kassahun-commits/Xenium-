#!/usr/bin/env python3
"""
DE within Neurons and Astrocytes: EtOH_MCT1i (test) vs EtOH_veh (reference),
using combined Slide A + Slide B data (n=3 vs n=2).

Reports:
  * Top genes UP in EtOH_MCT1i  (i.e. higher with drug than with vehicle)
  * Top genes UP in EtOH_veh    (i.e. higher with alcohol alone)
For each of: Neurons, Astrocytes.

Stats note: this is cell-level Wilcoxon. With 2 vs 3 biological replicates,
p-values are dominated by within-mouse cell variability. Direction and
magnitude (log2FC) are more reliable than raw p-values. For a publishable
analysis, run pseudobulk (per-sample) DE.
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

# (Same LABEL_MAP and MARKERS as the combined curated dotplot script.)
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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--slide-a-dir', required=True, type=Path)
    p.add_argument('--slide-a-ann', required=True, type=Path)
    p.add_argument('--slide-b-dir', required=True, type=Path)
    p.add_argument('--slide-b-ann', required=True, type=Path)
    p.add_argument('--out-dir', required=True, type=Path)
    p.add_argument('--label', default='SlidesAB')
    p.add_argument('--date', default=date.today().isoformat())
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(args.out_dir / f'{args.label}_DE_log.txt', mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== DE: EtOH_MCT1i vs EtOH_veh, Neurons + Astrocytes (combined slides) ===')

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

    # ---- DE: EtOH_MCT1i vs EtOH_veh per cell type ----
    for ct in ['Neuron', 'Astrocyte']:
        sub = adata[(adata.obs['celltype'].astype(str) == ct) &
                    (adata.obs['group'].astype(str).isin(['EtOH_veh', 'EtOH_MCT1i']))].copy()
        n_v = int((sub.obs['group'].astype(str) == 'EtOH_veh').sum())
        n_m = int((sub.obs['group'].astype(str) == 'EtOH_MCT1i').sum())
        s_v = sorted(sub.obs.loc[sub.obs['group'].astype(str) == 'EtOH_veh', 'sample_id'].unique().tolist())
        s_m = sorted(sub.obs.loc[sub.obs['group'].astype(str) == 'EtOH_MCT1i', 'sample_id'].unique().tolist())
        logging.info('%s: EtOH_veh=%d cells (samples %s), EtOH_MCT1i=%d cells (samples %s)',
                     ct, n_v, s_v, n_m, s_m)
        sub.obs['de_group'] = sub.obs['group'].astype(str)
        sc.tl.rank_genes_groups(sub, 'de_group', method='wilcoxon',
                                reference='EtOH_veh', n_genes=sub.shape[1])
        r = sub.uns['rank_genes_groups']
        de_df = pd.DataFrame({
            'gene': [str(g) for g in r['names']['EtOH_MCT1i']],
            'logfc': r['logfoldchanges']['EtOH_MCT1i'],
            'pval': r['pvals']['EtOH_MCT1i'],
            'padj': r['pvals_adj']['EtOH_MCT1i'],
            'score': r['scores']['EtOH_MCT1i'],
        })
        de_df['celltype'] = ct
        de_df['comparison'] = 'EtOH_MCT1i_vs_EtOH_veh'
        de_df['n_cells_EtOH_MCT1i'] = n_m
        de_df['n_cells_EtOH_veh'] = n_v
        de_df['n_samples_EtOH_MCT1i'] = len(s_m)
        de_df['n_samples_EtOH_veh'] = len(s_v)
        out_csv = args.out_dir / f'{args.label}_DE_{ct}_EtOH_MCT1i_vs_EtOH_veh_{args.date}.csv'
        de_df.to_csv(out_csv, index=False)
        logging.info('Wrote %s', out_csv)

        # Print top 25 each direction
        sig = de_df[(de_df['logfc'].abs() >= 1) & (de_df['padj'] < 1e-3)]
        up_in_mct = sig[sig['logfc'] > 0].sort_values('logfc', ascending=False).head(25)
        up_in_veh = sig[sig['logfc'] < 0].sort_values('logfc').head(25)
        logging.info('--- %s --- top 25 UP in EtOH_MCT1i (logfc>=1, padj<1e-3):\n%s',
                     ct, up_in_mct[['gene', 'logfc', 'padj']].to_string(index=False))
        logging.info('--- %s --- top 25 UP in EtOH_veh (logfc<=-1, padj<1e-3):\n%s',
                     ct, up_in_veh[['gene', 'logfc', 'padj']].to_string(index=False))

        # Volcano plot
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        x = de_df['logfc'].values
        y = -np.log10(np.clip(de_df['padj'].values, 1e-300, 1))
        s = (de_df['logfc'].abs() >= 1) & (de_df['padj'] < 1e-3)
        ax.scatter(x[~s], y[~s], s=6, c='lightgray', alpha=0.5, edgecolor='none')
        ax.scatter(x[s & (x > 0)], y[s & (x > 0)], s=14, c='#d62728', alpha=0.8,
                   edgecolor='none', label='UP in EtOH_MCT1i')
        ax.scatter(x[s & (x < 0)], y[s & (x < 0)], s=14, c='#1f77b4', alpha=0.8,
                   edgecolor='none', label='UP in EtOH_veh')
        top_up = de_df[s & (de_df['logfc'] > 0)].nlargest(15, 'logfc')
        top_dn = de_df[s & (de_df['logfc'] < 0)].nsmallest(15, 'logfc')
        for _, row in pd.concat([top_up, top_dn]).iterrows():
            ax.annotate(row['gene'], (row['logfc'], -np.log10(max(row['padj'], 1e-300))),
                        fontsize=7, alpha=0.85)
        ax.axhline(-np.log10(0.001), ls='--', color='k', lw=0.5, alpha=0.4)
        ax.axvline(1, ls='--', color='k', lw=0.5, alpha=0.4)
        ax.axvline(-1, ls='--', color='k', lw=0.5, alpha=0.4)
        ax.set_xlabel('log2 fold-change (EtOH_MCT1i / EtOH_veh)')
        ax.set_ylabel('-log10 adj p-value')
        ax.set_title(f'{ct} — EtOH_MCT1i (n={len(s_m)} mice) vs EtOH_veh (n={len(s_v)} mice)\n'
                     'Cell-level Wilcoxon (combined Slides A+B)')
        ax.legend(fontsize=8, frameon=False)
        plt.tight_layout()
        out_pdf = args.out_dir / f'{args.label}_DE_{ct}_EtOH_MCT1i_vs_EtOH_veh_{args.date}.pdf'
        plt.savefig(out_pdf, bbox_inches='tight')
        plt.close()
        logging.info('Wrote %s', out_pdf)

    logging.info('Done.')


if __name__ == '__main__':
    main()
