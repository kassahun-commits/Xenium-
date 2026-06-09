#!/usr/bin/env python3
"""
Cell-level Wilcoxon differential expression + volcano plots on the nucleus-only
combined object (Slides A+B), per cell type, each treatment group vs H2O_veh.

WHY THIS VERSION EXISTS / IMPORTANT CAVEAT
------------------------------------------
This reproduces the *exact* DE method used in the earlier Slide-B grant script
(`slideB_alcohol_prelim.py`): scanpy `rank_genes_groups(method='wilcoxon')` run
on individual cells, reference = H2O_veh, sig = padj < 0.01 & |log2FC| > 0.5.

It is provided for CONTINUITY with figures already used in a grant. It is
**EXPLORATORY ONLY** and is *not* the statistically preferred analysis:

  - Cell-level Wilcoxon treats every cell as an independent sample. With only
    n=2-3 mice per group, this is PSEUDOREPLICATION: p-values reflect the number
    of cells sequenced, not whether an effect reproduces across animals, so they
    are massively anti-conservative (inflated significance). See Squair et al.
    2021, Nat Commun, "Confronting false discoveries in single-cell DE".
  - The replicate-level pseudobulk DESeq2 version
    (`nucleus_pseudobulk_volcano_AB.py`) is the rigorous companion and should be
    the primary analysis. Use this Wilcoxon version only alongside it / for
    continuity.
  - Slide is also partially confounded with group (MAT2A_CM/OE are Slide-B-only;
    H2O_veh is split across slides) -- this affects BOTH methods and is worse,
    not better, under the cell-level test.

Same comparisons, thresholds, cell types, and output naming as the pseudobulk
script so the two can be placed side by side.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

CONTROL = 'H2O_veh'
GROUP_ORDER = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
               'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
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

# thresholds match the earlier grant script (slideB_alcohol_prelim.py)
LFC_THRESH = 0.5
PADJ_THRESH = 0.01
MIN_CELLS = 50   # skip a comparison with fewer than this many cells total


def wilcoxon_de(sub_ct, treatment, control=CONTROL):
    """Cell-level Wilcoxon (scanpy rank_genes_groups), treatment vs control.
    Mirrors slideB_alcohol_prelim.py exactly. Returns de_df or None."""
    sub = sub_ct[sub_ct.obs['group'].astype(str).isin([treatment, control])].copy()
    if sub.n_obs < MIN_CELLS:
        return None
    n_test = int((sub.obs['group'].astype(str) == treatment).sum())
    n_ref = int((sub.obs['group'].astype(str) == control).sum())
    if n_test < 1 or n_ref < 1:
        return None
    sub.obs['de_group'] = sub.obs['group'].astype(str)
    sc.tl.rank_genes_groups(sub, 'de_group', method='wilcoxon',
                            reference=control, n_genes=sub.n_vars)
    r = sub.uns['rank_genes_groups']
    de = pd.DataFrame({
        'gene': [str(g) for g in r['names'][treatment]],
        'log2FoldChange': r['logfoldchanges'][treatment],
        'pvalue': r['pvals'][treatment],
        'padj': r['pvals_adj'][treatment],
        'score': r['scores'][treatment],
    })
    de['comparison'] = f'{treatment}_vs_{control}'
    de['n_cells_test'] = n_test
    de['n_cells_ref'] = n_ref
    return de


def volcano_ax(ax, de, treatment, control, celltype):
    x = de['log2FoldChange'].values
    y = -np.log10(np.clip(de['padj'].values, 1e-300, 1))
    sig = (de['padj'] < PADJ_THRESH) & (de['log2FoldChange'].abs() > LFC_THRESH)
    up = sig & (de['log2FoldChange'] > 0)
    dn = sig & (de['log2FoldChange'] < 0)

    ax.scatter(x[~sig], y[~sig], s=6, c='#cccccc', alpha=0.5,
               edgecolor='none', rasterized=True)
    ax.scatter(x[up], y[up], s=10, c=GROUP_COLORS.get(treatment, '#d62728'),
               alpha=0.85, edgecolor='none', label=f'Up in {treatment}')
    ax.scatter(x[dn], y[dn], s=10, c=GROUP_COLORS.get(control, '#7f7f7f'),
               alpha=0.85, edgecolor='none', label=f'Up in {control}')

    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='black', lw=0.5, alpha=0.4)
    ax.axvline(LFC_THRESH, ls='--', color='black', lw=0.5, alpha=0.4)
    ax.axvline(-LFC_THRESH, ls='--', color='black', lw=0.5, alpha=0.4)

    top_up = de[up].nlargest(8, 'log2FoldChange')
    top_dn = de[dn].nsmallest(8, 'log2FoldChange')
    for _, rr in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(rr['gene'],
                    (rr['log2FoldChange'], -np.log10(max(rr['padj'], 1e-300))),
                    fontsize=5, ha='center', va='bottom')

    ax.set_title(f'{celltype}: {treatment} vs {control}\n'
                 f'EXPLORATORY cell-level Wilcoxon '
                 f'(n_up={int(up.sum())}, n_down={int(dn.sum())})', fontsize=8)
    ax.set_xlabel(f'log2 fold change ({treatment}/{control})', fontsize=7)
    ax.set_ylabel('-log10 adj. p', fontsize=7)
    ax.tick_params(labelsize=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--date', default='2026-06-09')
    ap.add_argument('--celltypes', nargs='+', default=CELLTYPES)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    D = args.date

    print(f'[load] {args.h5ad}')
    adata = ad.read_h5ad(args.h5ad)
    # rank_genes_groups runs on X (log-normalized) -- same as the grant script
    print(f'[note] DE on X (log-norm). Wilcoxon is cell-level -> EXPLORATORY ONLY.')

    all_results = []
    for ct in args.celltypes:
        sub_ct = adata[adata.obs['celltype'].astype(str) == ct].copy()
        print(f'\n[{ct}] {sub_ct.n_obs} cells')
        per_ct = {}
        for trt in TREATMENTS:
            de = wilcoxon_de(sub_ct, trt)
            if de is None:
                print(f'  [skip] {trt} vs {CONTROL}: too few cells')
                continue
            de['celltype'] = ct
            per_ct[trt] = de
            all_results.append(de)
            nsig = ((de['padj'] < PADJ_THRESH) &
                    (de['log2FoldChange'].abs() > LFC_THRESH)).sum()
            print(f'  [ok]  {trt} vs {CONTROL}: '
                  f'{de["n_cells_test"].iloc[0]} vs {de["n_cells_ref"].iloc[0]} cells; '
                  f'{nsig} sig (padj<{PADJ_THRESH}, |lfc|>{LFC_THRESH})')

            fig, ax = plt.subplots(figsize=(5, 4.5))
            volcano_ax(ax, de, trt, CONTROL, ct)
            ax.legend(fontsize=6, loc='upper right', frameon=False)
            plt.tight_layout()
            f1 = out_dir / f'Nucleus_WilcoxonDE_Volcano_{ct}_{trt}_vs_{CONTROL}_{D}.pdf'
            plt.savefig(f1, bbox_inches='tight', dpi=200)
            plt.close(fig)

        present = [t for t in TREATMENTS if t in per_ct]
        if present:
            ncol = 3
            nrow = int(np.ceil(len(present) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.3 * nrow),
                                     squeeze=False)
            for i, trt in enumerate(present):
                volcano_ax(axes[i // ncol][i % ncol], per_ct[trt], trt, CONTROL, ct)
            for j in range(len(present), nrow * ncol):
                axes[j // ncol][j % ncol].axis('off')
            fig.suptitle(f'{ct} — cell-level Wilcoxon vs {CONTROL} '
                         f'(|log2FC|>{LFC_THRESH}, padj<{PADJ_THRESH}) — EXPLORATORY ONLY',
                         fontsize=11)
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            f2 = out_dir / f'Nucleus_WilcoxonDE_Volcano_{ct}_AllGroups_{D}.pdf'
            plt.savefig(f2, bbox_inches='tight', dpi=200)
            plt.close(fig)
            print(f'  [fig] {f2.name}')

    if all_results:
        big = pd.concat(all_results, ignore_index=True)
        big['significant'] = (big['padj'] < PADJ_THRESH) & \
                             (big['log2FoldChange'].abs() > LFC_THRESH)
        cols = ['celltype', 'comparison', 'gene', 'log2FoldChange', 'score',
                'pvalue', 'padj', 'n_cells_test', 'n_cells_ref', 'significant']
        big = big[[c for c in cols if c in big.columns]]
        res_csv = out_dir / f'Nucleus_WilcoxonDE_results_{D}.csv'
        big.to_csv(res_csv, index=False)
        print(f'\n[csv] {res_csv.name}  ({len(big)} rows)')
    print('[done]')


if __name__ == '__main__':
    main()
