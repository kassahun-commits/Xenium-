#!/usr/bin/env python3
"""
Compare how the CELL-TYPING method changes a differential-expression result,
holding the test fixed.

Same contrast (EtOH_veh vs H2O_veh), same cell types (Neuron, Astrocyte), same
test (cell-level Wilcoxon, scanpy rank_genes_groups on log-normalised data),
same significance thresholds. The ONLY thing that differs between the two
columns of the figure is how each cell was assigned a type:

  LEFT  "marker-typed"  : our V1 DE pipeline — score canonical marker sets per
                          cell, label by the top-scoring set
                          (results read from PI_Master_NvsA_V1_DE CSV).
  RIGHT "leiden-typed"  : the 10x workshop pipeline — unsupervised leiden
                          clusters, then clusters grouped into Neuron / Astrocyte
                          by their marker genes (computed here from the processed
                          h5ad's lognorm layer).

Output: a 2x2 volcano grid (rows = cell type, cols = typing method) + the
underlying DE table as CSV. Editable-text PDF (fonttype 42). No hardcoded paths.
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

TEST, REF = 'EtOH_veh', 'H2O_veh'   # overridden by CLI in main()
HIGHLIGHT = None                     # optional gene to flag in every panel
CELLTYPES = ['Neuron', 'Astrocyte']
# leiden clusters -> cell type, from the 10x-run marker genes
# (Cnp/Mag/Mog=oligo, Gad1/Gad2=inhib neuron, Slc1a3/Gja1=astro, etc.)
LEIDEN_TO_CT = {
    'Neuron':    ['3', '5', '8', '10', '11', '12', '13', '14', '17', '18'],
    'Astrocyte': ['1', '15'],
}
# same thresholds as the Wilcoxon analysis elsewhere in the deck
LFC_THRESH = 0.5
PADJ_THRESH = 1e-3

CT_COLOR = {'Neuron': '#1f77b4', 'Astrocyte': '#2ca02c'}


def wilcoxon_de(adata_ct):
    """Cell-level Wilcoxon EtOH_veh vs H2O_veh on the lognorm layer."""
    sub = adata_ct[adata_ct.obs['group'].astype(str).isin([TEST, REF])].copy()
    sub.X = sub.layers['lognorm']
    sub.obs['de_group'] = sub.obs['group'].astype(str)
    n_test = int((sub.obs['de_group'] == TEST).sum())
    n_ref = int((sub.obs['de_group'] == REF).sum())
    sc.tl.rank_genes_groups(sub, 'de_group', method='wilcoxon',
                            reference=REF, n_genes=sub.n_vars)
    r = sub.uns['rank_genes_groups']
    de = pd.DataFrame({
        'gene': [str(g) for g in r['names'][TEST]],
        'logfc': r['logfoldchanges'][TEST],
        'padj': r['pvals_adj'][TEST],
    })
    return de, n_test, n_ref


def volcano(ax, de, title, accent):
    de = de.dropna(subset=['logfc', 'padj'])
    x = de['logfc'].values
    y = -np.log10(np.clip(de['padj'].values, 1e-300, 1))
    sig = (de['padj'] < PADJ_THRESH) & (de['logfc'].abs() > LFC_THRESH)
    up = sig & (de['logfc'] > 0)
    dn = sig & (de['logfc'] < 0)
    ax.scatter(x[~sig.values], y[~sig.values], s=6, c='#cccccc', alpha=0.45,
               edgecolor='none', rasterized=True)
    ax.scatter(x[up.values], y[up.values], s=11, c=accent, alpha=0.9,
               edgecolor='none')
    ax.scatter(x[dn.values], y[dn.values], s=11, c='#7f7f7f', alpha=0.9,
               edgecolor='none')
    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='k', lw=0.6, alpha=0.4)
    for v in (-LFC_THRESH, LFC_THRESH):
        ax.axvline(v, ls='--', color='k', lw=0.6, alpha=0.4)
    lab = pd.concat([de[up].nlargest(6, 'logfc'), de[dn].nsmallest(6, 'logfc')])
    for _, rr in lab.iterrows():
        if HIGHLIGHT and str(rr['gene']).lower() == HIGHLIGHT.lower():
            continue  # highlighted gene labelled separately below
        ax.annotate(rr['gene'],
                    (rr['logfc'], -np.log10(max(rr['padj'], 1e-300))),
                    fontsize=5.5, ha='center', va='bottom', color='#333')
    # flag the gene of interest (e.g. Grin2d) in every panel
    hl = None
    if HIGHLIGHT:
        hg = de[de['gene'].astype(str).str.lower() == HIGHLIGHT.lower()]
        if len(hg):
            rr = hg.iloc[0]
            hx = rr['logfc']
            hy = -np.log10(max(rr['padj'], 1e-300))
            hsig = bool((rr['padj'] < PADJ_THRESH) and
                        (abs(rr['logfc']) > LFC_THRESH))
            ax.scatter([hx], [hy], s=90, marker='*', c='#d62728',
                       edgecolor='k', linewidth=0.6, zorder=6)
            ax.annotate(f"{rr['gene']}  log2FC={rr['logfc']:.2f}, "
                        f"padj={rr['padj']:.1e}  [{'SIG' if hsig else 'n.s.'}]",
                        (hx, hy), fontsize=6.5, ha='left', va='bottom',
                        color='#d62728', fontweight='bold', zorder=7)
            hl = (float(rr['logfc']), float(rr['padj']), hsig)
    ax.set_title(f'{title}\n{int(up.sum())} up / {int(dn.sum())} down '
                 f'(padj<{PADJ_THRESH:g}, |log2FC|>{LFC_THRESH})', fontsize=9)
    ax.set_xlabel(f'log2 fold change ({TEST}/{REF})', fontsize=8)
    ax.set_ylabel('-log10 adj. p', fontsize=8)
    ax.tick_params(labelsize=7)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    return int(up.sum()), int(dn.sum()), set(de[sig]['gene']), hl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True, help='processed 10x-pipeline h5ad')
    ap.add_argument('--marker-de-csv', required=True,
                    help='PI_Master_NvsA_V1_DE Wilcoxon CSV (marker-typed)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--test', default='EtOH_veh', help='treatment group')
    ap.add_argument('--ref', default='H2O_veh', help='reference group')
    ap.add_argument('--highlight', default=None,
                    help='gene to flag in every panel (e.g. Grin2d)')
    ap.add_argument('--date', default='2026-06-09')
    args = ap.parse_args()

    global TEST, REF, HIGHLIGHT
    TEST, REF, HIGHLIGHT = args.test, args.ref, args.highlight

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    D = args.date

    # ---- marker-typed (already computed) ----
    mk = pd.read_csv(args.marker_de_csv)
    mk = mk[(mk['test'] == TEST) & (mk['reference'] == REF)].copy()
    mk = mk.rename(columns={'subtype': 'celltype'})

    # ---- leiden-typed (compute here) ----
    print(f'[load] {args.h5ad}')
    adata = ad.read_h5ad(args.h5ad)
    # free the big scaled X / counts; keep only lognorm + obs/var
    adata.X = adata.layers['lognorm']
    if 'counts' in adata.layers:
        del adata.layers['counts']
    leiden = adata.obs['leiden'].astype(str)
    ct_assign = pd.Series('Other', index=adata.obs_names)
    for ct, cls in LEIDEN_TO_CT.items():
        ct_assign[leiden.isin(cls).values] = ct
    adata.obs['celltype_leiden'] = ct_assign.values
    print('[leiden-typed] cell counts:',
          adata.obs['celltype_leiden'].value_counts().to_dict())

    rows = []
    leiden_de = {}
    for ct in CELLTYPES:
        sub = adata[adata.obs['celltype_leiden'] == ct]
        de, nt, nr = wilcoxon_de(sub)
        de['celltype'] = ct
        leiden_de[ct] = de
        leiden_de[ct + '_n'] = (nt, nr)
        print(f'  [leiden] {ct}: {nt} {TEST} vs {nr} {REF} cells')

    # ---- 2x2 volcano: rows = cell type, cols = typing method ----
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.0))
    summary = {}
    hl_report = {}
    for i, ct in enumerate(CELLTYPES):
        m = mk[mk['celltype'] == ct].rename(columns={'logfc': 'logfc'})
        nt_m = int(m['n_test_cells'].iloc[0]) if len(m) else 0
        nr_m = int(m['n_ref_cells'].iloc[0]) if len(m) else 0
        u1, d1, set_m, hl_m = volcano(
            axes[i][0], m[['gene', 'logfc', 'padj']],
            f'{ct} — marker-typed (n={nt_m} vs {nr_m})', CT_COLOR[ct])
        nt_l, nr_l = leiden_de[ct + '_n']
        u2, d2, set_l, hl_l = volcano(
            axes[i][1], leiden_de[ct][['gene', 'logfc', 'padj']],
            f'{ct} — leiden-typed (n={nt_l} vs {nr_l})', CT_COLOR[ct])
        inter = len(set_m & set_l)
        uni = len(set_m | set_l)
        jac = (inter / uni) if uni else float('nan')
        summary[ct] = (u1 + d1, u2 + d2, inter, jac)
        hl_report[ct] = (hl_m, hl_l)
        axes[i][0].text(-0.18, 0.5, ct, transform=axes[i][0].transAxes,
                        rotation=90, va='center', ha='center',
                        fontsize=13, fontweight='bold', color=CT_COLOR[ct])

    fig.suptitle(f'Same test (Wilcoxon), same contrast ({TEST} vs {REF}) — '
                 f'does the cell-typing method change the DE result?',
                 fontsize=12, fontweight='bold')
    cap = '   |   '.join(
        f'{ct}: marker {summary[ct][0]} sig vs leiden {summary[ct][1]} sig; '
        f'{summary[ct][2]} shared (Jaccard {summary[ct][3]:.2f})'
        for ct in CELLTYPES)
    fig.text(0.5, 0.005, cap, ha='center', fontsize=8.5, color='#333')
    plt.tight_layout(rect=[0.02, 0.02, 1, 0.96])
    pdf = out / f'V1_celltyping_compare_Wilcoxon_{TEST}_vs_{REF}_{D}.pdf'
    png = out / f'V1_celltyping_compare_Wilcoxon_{TEST}_vs_{REF}_{D}.png'
    fig.savefig(pdf, bbox_inches='tight')
    fig.savefig(png, dpi=180, bbox_inches='tight')
    plt.close(fig)

    # ---- source data ----
    mk_out = mk[['celltype', 'gene', 'logfc', 'padj']].copy()
    mk_out['method'] = 'marker-typed'
    li = pd.concat([leiden_de[ct] for ct in CELLTYPES], ignore_index=True)
    li = li[['celltype', 'gene', 'logfc', 'padj']].copy()
    li['method'] = 'leiden-typed'
    big = pd.concat([mk_out, li], ignore_index=True)
    big['significant'] = (big['padj'] < PADJ_THRESH) & \
                         (big['logfc'].abs() > LFC_THRESH)
    csv = out / f'V1_celltyping_compare_Wilcoxon_{TEST}_vs_{REF}_{D}.csv'
    big.to_csv(csv, index=False)

    print(f'[saved] {pdf.name}')
    print(f'[saved] {png.name}')
    print(f'[saved] {csv.name}  ({len(big)} rows)')
    for ct in CELLTYPES:
        print(f'  {ct}: marker {summary[ct][0]} sig, leiden {summary[ct][1]} '
              f'sig, {summary[ct][2]} shared, Jaccard {summary[ct][3]:.2f}')
    if HIGHLIGHT:
        print(f'\n[highlight] {HIGHLIGHT} ({TEST} vs {REF}):')
        for ct in CELLTYPES:
            for meth, hl in zip(('marker-typed', 'leiden-typed'), hl_report[ct]):
                if hl is None:
                    print(f'  {ct:>9} {meth:>13}: not in gene set')
                else:
                    lfc, padj, sig = hl
                    print(f'  {ct:>9} {meth:>13}: log2FC={lfc:+.3f}, '
                          f'padj={padj:.2e}  -> {"SIGNIFICANT" if sig else "n.s."}'
                          f'  (cut: padj<{PADJ_THRESH:g}, |log2FC|>{LFC_THRESH})')


if __name__ == '__main__':
    main()
