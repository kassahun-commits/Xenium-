#!/usr/bin/env python3
"""
Excitatory-neuron-only analysis: EtOH_MCT1i vs EtOH_veh.

Pipeline:
  1. Combine Slide A + Slide B, ROI assignment, label map (same as before).
  2. Cell-type into Neuron/Astrocyte/Oligo/Microglia with marker modules.
  3. Within "Neuron" cells, sub-classify EXCITATORY neurons using a panel
     of canonical excitatory markers (Slc17a7, Camk2a/b, Satb2, Tbr1, Neurod6).
     Threshold on the gene-score (> 0 = excitatory).
  4. Restrict to EtOH_veh and EtOH_MCT1i groups within excitatory neurons.
  5. Run cell-level Wilcoxon DE.
  6. Save a ranked DE table + a box-and-whisker grid of the top ~24 changed
     genes (matched-width boxes per condition, per-mouse mean overlay).

Outputs:
  - SlidesAB_ExcitatoryNeuron_EtOH_MCT1i_vs_EtOH_veh_DE_<date>.csv
  - SlidesAB_ExcitatoryNeuron_EtOH_MCT1i_vs_EtOH_veh_BoxPlots_<date>.pdf
  - SlidesAB_ExcitatoryNeuron_cellcounts_per_sample_<date>.csv
  - SlidesAB_ExcitatoryNeuron_log.txt
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

# Canonical excitatory markers (multi-marker -> robust)
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']

GROUPS_FOCUS = ['EtOH_veh', 'EtOH_MCT1i']


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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--n-top-genes', type=int, default=24,
                    help='How many top up + top down to plot')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_ExcitatoryNeuron_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Excitatory-neuron only: EtOH_MCT1i vs EtOH_veh ===')

    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    logging.info('Combined: %s cells x %s genes', f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Cell-type via marker modules
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

    # Score excitatory markers within all cells (we'll only use within neurons)
    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    logging.info('Excitatory markers in panel: %s', excit_present)
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory',
                      random_state=0, n_bins=25)

    # Excitatory subtype: cells with celltype==Neuron AND score_Excitatory > 0
    is_excit = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Excitatory'] > 0))
    adata.obs['neuron_subtype'] = np.where(
        adata.obs['celltype'].astype(str) == 'Neuron',
        np.where(is_excit, 'Excitatory', 'Other_Neuron'),
        adata.obs['celltype'].astype(str))
    logging.info('Cell subtype counts:\n%s',
                 adata.obs['neuron_subtype'].value_counts().to_string())

    # Per-sample excitatory cell counts
    excit_only = adata[is_excit].copy()
    logging.info('Excitatory neurons: %s cells', f'{excit_only.shape[0]:,}')
    samp_counts = (excit_only.obs.groupby(['group', 'sample_id'], observed=True)
                   .size().reset_index(name='n_excit'))
    samp_counts.to_csv(args.out_dir / f'{args.label}_ExcitatoryNeuron_cellcounts_per_sample_{args.date}.csv',
                       index=False)
    logging.info('Per-sample excitatory cell counts:\n%s', samp_counts.to_string(index=False))

    # Subset to EtOH_veh + EtOH_MCT1i
    sub = excit_only[excit_only.obs['group'].astype(str).isin(GROUPS_FOCUS)].copy()
    nv = int((sub.obs['group'].astype(str) == 'EtOH_veh').sum())
    nm = int((sub.obs['group'].astype(str) == 'EtOH_MCT1i').sum())
    s_v = sorted(sub.obs.loc[sub.obs['group'].astype(str) == 'EtOH_veh', 'sample_id'].unique().tolist())
    s_m = sorted(sub.obs.loc[sub.obs['group'].astype(str) == 'EtOH_MCT1i', 'sample_id'].unique().tolist())
    logging.info('Excitatory in EtOH_veh: %d cells (%s)', nv, s_v)
    logging.info('Excitatory in EtOH_MCT1i: %d cells (%s)', nm, s_m)

    # DE
    sub.obs['de_group'] = sub.obs['group'].astype(str)
    sc.tl.rank_genes_groups(sub, 'de_group', method='wilcoxon',
                            reference='EtOH_veh', n_genes=sub.shape[1])
    r = sub.uns['rank_genes_groups']
    de = pd.DataFrame({
        'gene': [str(g) for g in r['names']['EtOH_MCT1i']],
        'logfc': r['logfoldchanges']['EtOH_MCT1i'],
        'pval': r['pvals']['EtOH_MCT1i'],
        'padj': r['pvals_adj']['EtOH_MCT1i'],
        'score': r['scores']['EtOH_MCT1i'],
    })
    de['n_excit_EtOH_veh'] = nv
    de['n_excit_EtOH_MCT1i'] = nm
    de['n_mice_EtOH_veh'] = len(s_v)
    de['n_mice_EtOH_MCT1i'] = len(s_m)

    de_csv = args.out_dir / f'{args.label}_ExcitatoryNeuron_EtOH_MCT1i_vs_EtOH_veh_DE_{args.date}.csv'
    de.to_csv(de_csv, index=False)
    logging.info('Wrote %s', de_csv)

    # Pick top genes (significant + biggest |logFC|)
    sig = de[(de['logfc'].abs() >= 1) & (de['padj'] < 1e-3)].copy()
    top_up = sig[sig['logfc'] > 0].sort_values('logfc', ascending=False).head(args.n_top_genes)
    top_dn = sig[sig['logfc'] < 0].sort_values('logfc').head(args.n_top_genes)
    top_genes = top_up['gene'].tolist() + top_dn['gene'].tolist()
    logging.info('Top UP in MCT1i (%d): %s', len(top_up), top_up['gene'].tolist())
    logging.info('Top UP in EtOH_veh (%d): %s', len(top_dn), top_dn['gene'].tolist())

    # ---- Box plot grid ----
    np.random.seed(0)
    n_genes = len(top_genes)
    n_cols = 6
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(2.4 * n_cols, 2.0 * n_rows),
                              squeeze=False)

    # Per-cell long-form for plotting
    g_col = sub.obs['group'].astype(str).values
    s_col = sub.obs['sample_id'].astype(str).values

    for i, gene in enumerate(top_genes):
        r_, c_ = divmod(i, n_cols)
        ax = axes[r_, c_]
        if gene not in sub.var_names:
            ax.set_title(gene + ' (missing)', fontsize=8)
            ax.axis('off'); continue
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        d = pd.DataFrame({'expr': v, 'group': g_col, 'sample_id': s_col})

        all_vals = d['expr'].values
        q90 = float(np.percentile(all_vals, 90)) if len(all_vals) else 0.0
        use_strip = q90 < 0.5

        y_max = 0.0
        positions = {'EtOH_veh': 0, 'EtOH_MCT1i': 1}
        colors = {'EtOH_veh': '#ff7f0e', 'EtOH_MCT1i': '#2ca02c'}
        for grp in GROUPS_FOCUS:
            dg = d[d['group'] == grp]
            if len(dg) == 0: continue
            pos = positions[grp]
            sm = dg.groupby('sample_id')['expr'].mean()
            y_max = max(y_max, float(sm.max()) if len(sm) else 0.0)

            if not use_strip:
                ax.boxplot([dg['expr'].values], positions=[pos], widths=0.5,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color='black', lw=0.8),
                            whiskerprops=dict(color='black', lw=0.6),
                            capprops=dict(color='black', lw=0.6),
                            boxprops=dict(facecolor=colors[grp],
                                          edgecolor='black', linewidth=0.5, alpha=0.55))
                jitter = np.random.uniform(-0.08, 0.08, size=len(sm))
                ax.scatter([pos] * len(sm) + jitter, sm.values,
                            c='black', s=20, zorder=5,
                            edgecolor='white', linewidth=0.5)
                y_max = max(y_max, float(np.percentile(dg['expr'].values, 95)))
            else:
                jitter = np.random.uniform(-0.12, 0.12, size=len(sm))
                ax.scatter([pos] * len(sm) + jitter, sm.values,
                            c=colors[grp], s=30, zorder=4,
                            edgecolor='black', linewidth=0.5, alpha=0.85)

        # Annotate logFC + padj
        gene_de = de[de['gene'] == gene]
        if len(gene_de):
            lfc = float(gene_de['logfc'].iloc[0])
            padj = float(gene_de['padj'].iloc[0])
            ann = f'lfc={lfc:+.1f}\npadj={padj:.1e}'
            ax.text(0.98, 0.97, ann, transform=ax.transAxes,
                    fontsize=6, va='top', ha='right',
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='gray',
                              alpha=0.85))
        ax.set_xticks([0, 1])
        ax.set_xticklabels(GROUPS_FOCUS, rotation=20, ha='right', fontsize=7)
        ax.set_title(gene, fontsize=10, fontweight='bold')
        ax.tick_params(axis='y', labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylim(-0.05 * max(y_max, 0.05), max(y_max * 1.18, 0.1))
        if use_strip:
            ax.text(0.02, 0.97, 'low expr',
                    transform=ax.transAxes, fontsize=6, va='top', ha='left',
                    style='italic', color='#777')

    for j in range(n_genes, n_rows * n_cols):
        r_, c_ = divmod(j, n_cols)
        axes[r_, c_].axis('off')

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor='#ff7f0e', edgecolor='black', alpha=0.55, label='EtOH_veh'),
        Patch(facecolor='#2ca02c', edgecolor='black', alpha=0.55, label='EtOH_MCT1i'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                   markersize=6, label='per-mouse mean'),
    ]
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=3, frameon=False, fontsize=10)
    fig.suptitle(
        f'{args.label}: Excitatory neurons only — EtOH_MCT1i vs EtOH_veh\n'
        f'top {len(top_up)} UP (in MCT1i) + top {len(top_dn)} UP (in EtOH_veh) genes by |log2FC|; '
        f'EtOH_veh n={len(s_v)} mice / {nv:,} cells; EtOH_MCT1i n={len(s_m)} mice / {nm:,} cells',
        fontsize=11, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.99])

    pdf_out = args.out_dir / f'{args.label}_ExcitatoryNeuron_EtOH_MCT1i_vs_EtOH_veh_BoxPlots_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', pdf_out)


if __name__ == '__main__':
    main()
