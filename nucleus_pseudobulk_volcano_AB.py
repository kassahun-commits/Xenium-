#!/usr/bin/env python3
"""
Pseudobulk DESeq2 differential expression + volcano plots on the nucleus-only
combined object (Slides A+B), per cell type, each treatment group vs H2O_veh.

WHY PSEUDOBULK (not cell-level Wilcoxon):
  The earlier Slide-B-only script (slideB_alcohol_prelim.py) used cell-level
  Wilcoxon because, on Slide B alone, each alcohol group was a single tissue
  punch (n=1 biological replicate) -- there was no replicate structure to model,
  so that test was explicitly labelled "EXPLORATORY ONLY" in its code.
  Now that both slides are combined we have n=2-3 mice per group, so the
  statistically correct approach is pseudobulk: sum raw nuclear counts within
  each (cell type, biological replicate), then run DESeq2 across replicates.
  The mouse -- not the cell -- is the experimental unit.

DESIGN NOTE / LIMITATION (documented in outputs):
  - Per-group n is small (2-3). Power is limited; absence of significant genes
    is not evidence of no effect.
  - Slide is partially confounded with group: MAT2A_CM / MAT2A_OE replicates are
    all on Slide B, while H2O_veh is split (1 rep per slide). With only 2-3
    replicates per group there are not enough degrees of freedom to safely add
    slide as a covariate, so a simple ~group design is used. Interpret MAT2A
    contrasts with extra caution (group effect and slide effect are not
    separable for those).

Inputs read from CLI (no hardcoded paths). Outputs: per-comparison volcano PDFs
(editable text, fonttype 42), a combined multi-panel PDF per cell type, and a
source-data CSV of full DESeq2 results for every comparison.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ---- project constants (match existing pipeline) -------------------------
CONTROL = 'H2O_veh'
GROUP_ORDER = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
               'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
# all non-control groups, compared vs H2O_veh
TREATMENTS = [g for g in GROUP_ORDER if g != CONTROL]
CELLTYPES = ['Neuron', 'Astrocyte']

GROUP_COLORS = {
    'H2O_veh':     '#7f7f7f',
    'H2O_MCT1i':   '#17becf',
    'EtOH_veh':    '#ff7f0e',
    'EtOH_MCT1i':  '#2ca02c',
    'ChronicEtOH': '#d62728',
    'MAT2A_CM':    '#9467bd',
    'MAT2A_OE':    '#8c564b',
}

# DE thresholds
LFC_THRESH = 0.5          # |log2 fold change| cutoff (user requested 0.5/-0.5)
PADJ_THRESH = 0.05        # adjusted p-value cutoff (standard for pseudobulk DESeq2)
MIN_TOTAL_COUNT = 10      # gene-level pre-filter: total counts across samples in a comparison
MIN_CELLS_PER_SAMPLE = 25 # drop a pseudobulk replicate if it has fewer cells of this type


# ---- pseudobulk construction --------------------------------------------
def make_pseudobulk(adata_ct):
    """Sum raw counts (layers['counts']) per biological replicate (sample_id).
    Returns (counts_df [samples x genes, int], meta_df [sample_id, group, n_cells])."""
    counts = adata_ct.layers['counts']
    samples = adata_ct.obs['sample_id'].astype(str).values
    uniq = sorted(pd.unique(samples))
    genes = list(adata_ct.var_names)
    mat = np.zeros((len(uniq), adata_ct.n_vars), dtype=np.int64)
    n_cells = []
    groups = []
    for i, s in enumerate(uniq):
        m = samples == s
        sub = counts[m]
        sub = np.asarray(sub.todense()) if hasattr(sub, 'todense') else np.asarray(sub)
        mat[i] = np.rint(sub.sum(axis=0)).astype(np.int64)
        n_cells.append(int(m.sum()))
        groups.append(adata_ct.obs.loc[m, 'group'].astype(str).iloc[0])
    counts_df = pd.DataFrame(mat, index=uniq, columns=genes)
    meta_df = pd.DataFrame({'sample_id': uniq, 'group': groups, 'n_cells': n_cells},
                           index=uniq)
    return counts_df, meta_df


def run_deseq(counts_df, meta_df, treatment, control=CONTROL):
    """Run DESeq2 for treatment vs control on the two relevant groups.
    Returns results_df (index=gene) or None if a group is missing/too small."""
    keep_samples = meta_df.index[meta_df['group'].isin([treatment, control])]
    meta = meta_df.loc[keep_samples].copy()
    # require at least 2 replicates per side for a meaningful test
    g_counts = meta['group'].value_counts()
    if g_counts.get(treatment, 0) < 2 or g_counts.get(control, 0) < 2:
        return None, meta
    counts = counts_df.loc[keep_samples]
    # gene pre-filter: total counts across these samples
    counts = counts.loc[:, counts.sum(axis=0) >= MIN_TOTAL_COUNT]
    # pydeseq2 replaces underscores in factor levels with hyphens internally;
    # sanitize up front so the contrast levels match.
    ctrl_s, trt_s = control.replace('_', '-'), treatment.replace('_', '-')
    meta['group'] = pd.Categorical([g.replace('_', '-') for g in meta['group']],
                                   categories=[ctrl_s, trt_s])

    dds = DeseqDataSet(
        counts=counts,
        metadata=meta,
        design_factors='group',
        ref_level=['group', ctrl_s],
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['group', trt_s, ctrl_s], quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    res['gene'] = res.index
    res['celltype'] = None
    res['comparison'] = f'{treatment}_vs_{control}'
    return res, meta


# ---- plotting ------------------------------------------------------------
def volcano_ax(ax, res, treatment, control, celltype):
    """Draw one volcano on ax. res must have log2FoldChange + padj."""
    d = res.dropna(subset=['log2FoldChange', 'padj']).copy()
    # cap padj=0 to avoid inf on -log10
    minp = d.loc[d['padj'] > 0, 'padj'].min() if (d['padj'] > 0).any() else 1e-300
    d['padj_plot'] = d['padj'].clip(lower=minp)
    d['neglog10'] = -np.log10(d['padj_plot'])
    sig = (d['padj'] < PADJ_THRESH) & (d['log2FoldChange'].abs() > LFC_THRESH)
    up = sig & (d['log2FoldChange'] > 0)
    dn = sig & (d['log2FoldChange'] < 0)
    ns = ~sig

    ax.scatter(d.loc[ns, 'log2FoldChange'], d.loc[ns, 'neglog10'],
               s=6, c='#cccccc', alpha=0.6, edgecolor='none', rasterized=True)
    ax.scatter(d.loc[up, 'log2FoldChange'], d.loc[up, 'neglog10'],
               s=10, c=GROUP_COLORS.get(treatment, '#d62728'), alpha=0.85,
               edgecolor='none', label=f'Up in {treatment}')
    ax.scatter(d.loc[dn, 'log2FoldChange'], d.loc[dn, 'neglog10'],
               s=10, c=GROUP_COLORS.get(control, '#7f7f7f'), alpha=0.85,
               edgecolor='none', label=f'Up in {control}')

    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='black', lw=0.5, alpha=0.4)
    ax.axvline(LFC_THRESH, ls='--', color='black', lw=0.5, alpha=0.4)
    ax.axvline(-LFC_THRESH, ls='--', color='black', lw=0.5, alpha=0.4)

    # label top genes by padj among significant
    top = d.loc[sig].nsmallest(8, 'padj')
    for _, r in top.iterrows():
        ax.annotate(r['gene'], (r['log2FoldChange'], r['neglog10']),
                    fontsize=5, ha='center', va='bottom')

    ax.set_title(f'{celltype}: {treatment} vs {control}\n'
                 f'(n_up={int(up.sum())}, n_down={int(dn.sum())})', fontsize=8)
    ax.set_xlabel('log2 fold change', fontsize=7)
    ax.set_ylabel('-log10 adj. p', fontsize=7)
    ax.tick_params(labelsize=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True, help='combined nucleus-only h5ad')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--date', default='2026-06-09')
    ap.add_argument('--celltypes', nargs='+', default=CELLTYPES)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    D = args.date

    print(f'[load] {args.h5ad}')
    adata = ad.read_h5ad(args.h5ad)
    assert 'counts' in adata.layers, 'expected raw counts in layers["counts"]'

    all_results = []
    pb_meta_records = []

    for ct in args.celltypes:
        sub = adata[adata.obs['celltype'] == ct].copy()
        print(f'\n[{ct}] {sub.n_obs} cells')
        counts_df, meta_df = make_pseudobulk(sub)
        # record / drop tiny pseudobulk replicates
        for _, r in meta_df.iterrows():
            pb_meta_records.append({'celltype': ct, **r.to_dict()})
        small = meta_df.index[meta_df['n_cells'] < MIN_CELLS_PER_SAMPLE]
        if len(small):
            print(f'  [warn] dropping {len(small)} replicate(s) with <{MIN_CELLS_PER_SAMPLE} '
                  f'{ct} cells: {list(small)}')
            counts_df = counts_df.drop(index=small)
            meta_df = meta_df.drop(index=small)

        per_ct_res = {}
        ng_all = meta_df['group'].value_counts()
        for trt in TREATMENTS:
            res, meta = run_deseq(counts_df, meta_df, trt)
            if res is None:
                print(f'  [skip] {trt} vs {CONTROL}: insufficient replicates '
                      f'({trt}={ng_all.get(trt,0)}, {CONTROL}={ng_all.get(CONTROL,0)})')
                continue
            res['celltype'] = ct
            per_ct_res[trt] = res
            all_results.append(res)
            nsig = ((res['padj'] < PADJ_THRESH) &
                    (res['log2FoldChange'].abs() > LFC_THRESH)).sum()
            print(f'  [ok]  {trt} vs {CONTROL}: '
                  f'{trt}={ng_all.get(trt)}, {CONTROL}={ng_all.get(CONTROL)} reps; '
                  f'{nsig} sig genes (padj<{PADJ_THRESH}, |lfc|>{LFC_THRESH})')

            # single-comparison volcano
            fig, ax = plt.subplots(figsize=(5, 4.5))
            volcano_ax(ax, res, trt, CONTROL, ct)
            ax.legend(fontsize=6, loc='upper right', frameon=False)
            plt.tight_layout()
            f1 = out_dir / f'Nucleus_PseudobulkDE_Volcano_{ct}_{trt}_vs_{CONTROL}_{D}.pdf'
            plt.savefig(f1, bbox_inches='tight', dpi=200)
            plt.close(fig)

        # combined multi-panel per cell type
        present = [t for t in TREATMENTS if t in per_ct_res]
        if present:
            ncol = 3
            nrow = int(np.ceil(len(present) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.3 * nrow),
                                     squeeze=False)
            for i, trt in enumerate(present):
                volcano_ax(axes[i // ncol][i % ncol], per_ct_res[trt], trt, CONTROL, ct)
            for j in range(len(present), nrow * ncol):
                axes[j // ncol][j % ncol].axis('off')
            fig.suptitle(f'{ct} — pseudobulk DESeq2 vs {CONTROL} '
                         f'(|log2FC|>{LFC_THRESH}, padj<{PADJ_THRESH})', fontsize=11)
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            f2 = out_dir / f'Nucleus_PseudobulkDE_Volcano_{ct}_AllGroups_{D}.pdf'
            plt.savefig(f2, bbox_inches='tight', dpi=200)
            plt.close(fig)
            print(f'  [fig] {f2.name}')

    # ---- source data CSVs ----
    if all_results:
        big = pd.concat(all_results, ignore_index=True)
        cols = ['celltype', 'comparison', 'gene', 'baseMean', 'log2FoldChange',
                'lfcSE', 'stat', 'pvalue', 'padj']
        cols = [c for c in cols if c in big.columns]
        big = big[cols]
        big['significant'] = (big['padj'] < PADJ_THRESH) & \
                             (big['log2FoldChange'].abs() > LFC_THRESH)
        res_csv = out_dir / f'Nucleus_PseudobulkDE_results_{D}.csv'
        big.to_csv(res_csv, index=False)
        print(f'\n[csv] {res_csv.name}  ({len(big)} rows)')

    pb_df = pd.DataFrame(pb_meta_records)
    pb_csv = out_dir / f'Nucleus_Pseudobulk_replicates_per_celltype_{D}.csv'
    pb_df.to_csv(pb_csv, index=False)
    print(f'[csv] {pb_csv.name}  ({len(pb_df)} rows)')
    print('\n[done]')


if __name__ == '__main__':
    main()
