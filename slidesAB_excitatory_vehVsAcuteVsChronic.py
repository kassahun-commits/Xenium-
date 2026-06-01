#!/usr/bin/env python3
"""
Excitatory neurons only: H2O_veh vs EtOH_veh (acute) vs ChronicEtOH.

Identifies genes that change in either acute or chronic alcohol (vs control),
plots them as box-and-whisker grid with 3 boxes per panel (control / acute /
chronic).

Force-included regardless of significance: MCT family (Slc16a1=MCT1,
Slc16a7=MCT2, Slc16a3=MCT4). Gene labels use 'Slc... (alias)' for
well-known protein names.

Same plumbing as slidesAB_excitatory_MCT1i_vs_veh.py — combines Slide A + B,
ROI-assigns cells, marker-scored cell types, sub-classifies excitatory.
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
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']

GROUPS = ['H2O_veh', 'EtOH_veh', 'EtOH_MCT1i', 'ChronicEtOH']
GROUP_COLORS = {
    'H2O_veh':    '#888888',   # gray = control
    'EtOH_veh':   '#ff7f0e',   # orange = acute alcohol (veh)
    'EtOH_MCT1i': '#2ca02c',   # green = acute alcohol + MCT1i
    'ChronicEtOH':'#a83232',   # dark red = chronic
}

# Protein-name aliases for nicer labels
ALIAS = {
    'Slc16a1': 'MCT1',
    'Slc16a7': 'MCT2',
    'Slc16a3': 'MCT4',
    'Slc17a7': 'VGLUT1',
    'Slc17a6': 'VGLUT2',
    'Slc2a1':  'GLUT1',
    'Slc2a2':  'GLUT2',
    'Slc2a3':  'GLUT3',
    'Slc2a4':  'GLUT4',
    'Slc2a5':  'GLUT5',
    'Slc1a3':  'GLAST',
    'Slc1a1':  'EAAC1',
    'Slc1a2':  'GLT-1',
    'Slc32a1': 'VGAT',
    'Rbfox3':  'NeuN',
    'Crebbp':  'CBP',
    'Ep300':   'p300',
    'Kat2a':   'GCN5',
    'Kat2b':   'PCAF',
    'Kat5':    'Tip60',
}

# Force-include — always plot these even if not significant
FORCE_INCLUDE = ['Slc16a1', 'Slc16a7', 'Slc16a3']


def disp(gene: str) -> str:
    return f'{gene} ({ALIAS[gene]})' if gene in ALIAS else gene


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
                    help='How many top up + top down genes to plot per contrast')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_ExcitatoryNeuron_vehAcuteChronic_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Excitatory neurons: H2O_veh vs EtOH_veh vs ChronicEtOH ===')

    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    logging.info('Combined: %s cells x %s genes', f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Cell-type scoring
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

    # Excitatory subtype
    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory',
                      random_state=0, n_bins=25)
    is_excit = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Excitatory'] > 0))
    excit = adata[is_excit].copy()
    logging.info('Excitatory neurons total: %d', excit.shape[0])

    # Restrict to 3 groups of interest
    sub = excit[excit.obs['group'].astype(str).isin(GROUPS)].copy()
    sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str),
                                      categories=GROUPS, ordered=True)

    # Per-sample / per-group counts
    counts = (sub.obs.groupby(['group', 'sample_id'], observed=True)
              .size().reset_index(name='n_excit'))
    counts.to_csv(args.out_dir / f'{args.label}_Excit_vehAcuteChronic_cellcounts_{args.date}.csv',
                  index=False)
    logging.info('Per-sample excitatory counts:\n%s', counts.to_string(index=False))

    n_per_group = sub.obs['group'].value_counts()
    n_mice = {g: sub.obs.loc[sub.obs['group'] == g, 'sample_id'].nunique() for g in GROUPS}
    logging.info('Cells per group: %s', n_per_group.to_dict())
    logging.info('Mice per group: %s', n_mice)

    # Three pairwise DE vs H2O_veh: EtOH_veh, EtOH_MCT1i, ChronicEtOH
    sub.obs['de_group'] = sub.obs['group'].astype(str)
    de_combined = []
    for test in ['EtOH_veh', 'EtOH_MCT1i', 'ChronicEtOH']:
        d2 = sub[sub.obs['de_group'].isin([test, 'H2O_veh'])].copy()
        sc.tl.rank_genes_groups(d2, 'de_group', method='wilcoxon',
                                reference='H2O_veh', n_genes=d2.shape[1])
        r = d2.uns['rank_genes_groups']
        df = pd.DataFrame({
            'gene':  [str(g) for g in r['names'][test]],
            'logfc': r['logfoldchanges'][test],
            'pval':  r['pvals'][test],
            'padj':  r['pvals_adj'][test],
            'score': r['scores'][test],
        })
        df['contrast'] = f'{test}_vs_H2O_veh'
        de_combined.append(df)
    de_long = pd.concat(de_combined, ignore_index=True)
    de_long.to_csv(args.out_dir / f'{args.label}_Excit_vehAcuteChronic_DE_{args.date}.csv', index=False)
    logging.info('Wrote DE table.')

    # Pivot to wide for picking top hits
    wide = (de_long.pivot(index='gene', columns='contrast',
                          values=['logfc', 'padj']))
    # max |logfc| across the three contrasts; require sig in at least one
    lfc_acute = wide[('logfc', 'EtOH_veh_vs_H2O_veh')]
    lfc_mct   = wide[('logfc', 'EtOH_MCT1i_vs_H2O_veh')]
    lfc_chr   = wide[('logfc', 'ChronicEtOH_vs_H2O_veh')]
    padj_a    = wide[('padj',  'EtOH_veh_vs_H2O_veh')]
    padj_m    = wide[('padj',  'EtOH_MCT1i_vs_H2O_veh')]
    padj_c    = wide[('padj',  'ChronicEtOH_vs_H2O_veh')]
    max_abs   = pd.concat([lfc_acute.abs(), lfc_mct.abs(), lfc_chr.abs()], axis=1).max(axis=1)
    any_sig   = (padj_a < 1e-3) | (padj_m < 1e-3) | (padj_c < 1e-3)
    any_big   = (lfc_acute.abs() >= 1.0) | (lfc_mct.abs() >= 1.0) | (lfc_chr.abs() >= 1.0)
    hit = any_sig & any_big

    # signed: take the lfc with largest |value| among the three contrasts
    lfc_stack = pd.concat([lfc_acute.fillna(0), lfc_mct.fillna(0), lfc_chr.fillna(0)], axis=1)
    abs_stack = lfc_stack.abs()
    signed = lfc_stack.values[np.arange(len(lfc_stack)), abs_stack.values.argmax(axis=1)]
    signed = pd.Series(signed, index=wide.index)

    sig_df = wide.loc[hit].copy()
    sig_df['max_abs_lfc'] = max_abs[hit]
    sig_df['signed_lfc'] = signed[hit]
    top_up = sig_df[sig_df['signed_lfc'] > 0].sort_values('max_abs_lfc', ascending=False).head(args.n_top_genes)
    top_dn = sig_df[sig_df['signed_lfc'] < 0].sort_values('max_abs_lfc', ascending=False).head(args.n_top_genes)
    top_genes = top_up.index.tolist() + top_dn.index.tolist()
    logging.info('Top UP genes (acute or chronic vs H2O): %s', top_up.index.tolist())
    logging.info('Top DOWN genes (acute or chronic vs H2O): %s', top_dn.index.tolist())

    # Force-include MCTs
    for g in FORCE_INCLUDE:
        if g in adata.var_names and g not in top_genes:
            top_genes.append(g)
            logging.info('Force-included: %s', g)

    # ---- Box plot grid ----
    np.random.seed(0)
    n_genes = len(top_genes)
    n_cols = 6
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(2.6 * n_cols, 2.1 * n_rows),
                              squeeze=False)

    g_col = sub.obs['group'].astype(str).values
    s_col = sub.obs['sample_id'].astype(str).values

    for i, gene in enumerate(top_genes):
        r_, c_ = divmod(i, n_cols)
        ax = axes[r_, c_]
        if gene not in sub.var_names:
            ax.set_title(disp(gene) + ' (missing)', fontsize=8); ax.axis('off'); continue
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        d = pd.DataFrame({'expr': v, 'group': g_col, 'sample_id': s_col})

        all_vals = d['expr'].values
        q90 = float(np.percentile(all_vals, 90)) if len(all_vals) else 0.0
        use_strip = q90 < 0.5

        y_max = 0.0
        positions = {g: i for i, g in enumerate(GROUPS)}
        for grp in GROUPS:
            dg = d[d['group'] == grp]
            if len(dg) == 0: continue
            pos = positions[grp]
            sm = dg.groupby('sample_id')['expr'].mean()
            y_max = max(y_max, float(sm.max()) if len(sm) else 0.0)
            if not use_strip:
                ax.boxplot([dg['expr'].values], positions=[pos], widths=0.55,
                           patch_artist=True, showfliers=False,
                           medianprops=dict(color='black', lw=0.8),
                           whiskerprops=dict(color='black', lw=0.6),
                           capprops=dict(color='black', lw=0.6),
                           boxprops=dict(facecolor=GROUP_COLORS[grp],
                                         edgecolor='black', linewidth=0.5, alpha=0.6))
                jitter = np.random.uniform(-0.08, 0.08, size=len(sm))
                ax.scatter([pos] * len(sm) + jitter, sm.values,
                           c='black', s=18, zorder=5,
                           edgecolor='white', linewidth=0.5)
                y_max = max(y_max, float(np.percentile(dg['expr'].values, 95)))
            else:
                jitter = np.random.uniform(-0.12, 0.12, size=len(sm))
                ax.scatter([pos] * len(sm) + jitter, sm.values,
                           c=GROUP_COLORS[grp], s=28, zorder=4,
                           edgecolor='black', linewidth=0.5, alpha=0.85)

        # Annotate logFCs
        gene_de_a = de_long[(de_long['gene'] == gene) & (de_long['contrast'] == 'EtOH_veh_vs_H2O_veh')]
        gene_de_m = de_long[(de_long['gene'] == gene) & (de_long['contrast'] == 'EtOH_MCT1i_vs_H2O_veh')]
        gene_de_c = de_long[(de_long['gene'] == gene) & (de_long['contrast'] == 'ChronicEtOH_vs_H2O_veh')]
        ann_lines = []
        if len(gene_de_a):
            ann_lines.append(f"acute lfc={gene_de_a['logfc'].iloc[0]:+.1f}")
        if len(gene_de_m):
            ann_lines.append(f"+MCT1i lfc={gene_de_m['logfc'].iloc[0]:+.1f}")
        if len(gene_de_c):
            ann_lines.append(f"chronic lfc={gene_de_c['logfc'].iloc[0]:+.1f}")
        ax.text(0.98, 0.97, '\n'.join(ann_lines), transform=ax.transAxes,
                fontsize=6, va='top', ha='right',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='gray', alpha=0.85))

        ax.set_xticks(list(range(len(GROUPS))))
        ax.set_xticklabels(GROUPS, rotation=25, ha='right', fontsize=7)
        ax.set_title(disp(gene), fontsize=10, fontweight='bold')
        ax.tick_params(axis='y', labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.6, len(GROUPS) - 0.4)
        ax.set_ylim(-0.05 * max(y_max, 0.05), max(y_max * 1.18, 0.1))
        if use_strip:
            ax.text(0.02, 0.97, 'low expr',
                    transform=ax.transAxes, fontsize=6, va='top', ha='left',
                    style='italic', color='#777')
        if gene in FORCE_INCLUDE:
            ax.text(0.02, 0.02, 'force-included',
                    transform=ax.transAxes, fontsize=6, va='bottom', ha='left',
                    style='italic', color='#444')

    for j in range(n_genes, n_rows * n_cols):
        r_, c_ = divmod(j, n_cols)
        axes[r_, c_].axis('off')

    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor='black',
                            alpha=0.6, label=g) for g in GROUPS]
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                      markerfacecolor='black', markersize=6,
                                      label='per-mouse mean'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=4, frameon=False, fontsize=10)

    n_str = ', '.join([f"{g} n={n_mice[g]}" for g in GROUPS])
    fig.suptitle(
        f'{args.label}: Excitatory neurons — H2O_veh vs EtOH_veh vs EtOH_MCT1i vs ChronicEtOH\n'
        f'top {len(top_up)} UP + top {len(top_dn)} DOWN by max(|log2FC|) across the three contrasts vs H2O_veh; '
        f'plus MCT family (force-included). {n_str}',
        fontsize=11, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.99])

    pdf_out = args.out_dir / f'{args.label}_Excit_vehAcuteChronic_BoxPlots_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', pdf_out)


if __name__ == '__main__':
    main()
