#!/usr/bin/env python3
"""
Focused dot plot — same curated gene list, but filtered to genes that
visually CHANGE in either of:
  (a) Neuron vs Astrocyte (combined across groups), OR
  (b) EtOH_veh vs EtOH_MCT1i (within Neurons or within Astrocytes).

Uses the same combined Slide A + Slide B pipeline as
`slidesAB_combined_curated_dotplot.py`. Same layout (Neuron | Astrocyte
blocks, 7 condition columns each, categories on the right).
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
    # Slide B
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
    # Slide A
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
GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
          'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
CELLTYPES = ['Neuron', 'Astrocyte']

GENE_CATEGORIES = {
    'Cell-type identity': ['Rbfox3', 'Snap25', 'Syn1', 'Map2',
                           'Gfap', 'Aqp4', 'Slc1a3', 'S100b', 'Aldh1l1'],
    'Alcohol & acetaldehyde metab': ['Adh1', 'Adh5', 'Aldh2', 'Aldh1a1',
                                     'Cyp2e1', 'Cat', 'Akr1a1'],
    'Acetate / Ac-CoA metab': ['Acss1', 'Acss2', 'Acly', 'Acaca', 'Acacb'],
    'Lactate / MCT axis': ['Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb'],
    'Glucose transport': ['Slc2a1', 'Slc2a3', 'Slc2a4'],
    'Glycolysis / TCA': ['Gapdh', 'Pgk1', 'Pkm', 'Hk1', 'Hk2', 'Pfkm', 'Aldoa',
                         'Eno1', 'Eno2', 'Idh1', 'Idh2', 'Cs'],
    'SAM / MAT2A axis': ['Mat1a', 'Mat2a', 'Mat2b', 'Ahcy', 'Mthfr', 'Cbs'],
    'Glutathione / oxidative stress': ['Gclc', 'Gclm', 'Gpx1', 'Gpx4', 'Gss',
                                       'Gstm1', 'Sod1', 'Sod2', 'Hmox1', 'Nfe2l2'],
    'Stress / UPR': ['Atf3', 'Atf4', 'Ddit3', 'Hspa1a', 'Hsp90aa1', 'Hspb1', 'Trp53'],
    'Activity / IEGs': ['Fos', 'Jun', 'Junb', 'Arc', 'Egr1', 'Egr2',
                        'Npas4', 'Homer1', 'Bdnf'],
    'Glutamate signaling': ['Slc17a7', 'Slc17a6', 'Slc1a2', 'Grin1', 'Grin2a',
                            'Grin2b', 'Gria1', 'Gria2', 'Grm1', 'Grm2', 'Grm5',
                            'Camk2a', 'Camk2b'],
    'GABA signaling': ['Gad1', 'Gad2', 'Gabra1', 'Gabra2', 'Gabbr1', 'Slc32a1'],
    'Endocannabinoid': ['Cnr1', 'Faah', 'Mgll', 'Daglb'],
    'Neuroinflammation': ['Tnf', 'Il1b', 'Il6', 'Il10', 'Tlr4', 'Nfkb1', 'Stat3'],
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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--celltype-delta', type=float, default=0.5,
                    help='Min |mean_Neuron - mean_Astrocyte| in any group (log-norm)')
    ap.add_argument('--treatment-delta', type=float, default=0.3,
                    help='Min |mean_EtOH_veh - mean_EtOH_MCT1i| within a celltype (log-norm)')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout)])
    logging.info('=== Focused dot plot — combined slides ===')
    logging.info('Thresholds: celltype Δ=%.2f, treatment Δ=%.2f', args.celltype_delta, args.treatment_delta)

    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    logging.info('Combined: %s cells x %s genes', f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

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

    sub = adata[adata.obs['celltype'].astype(str).isin(CELLTYPES)].copy()
    sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str), categories=GROUPS, ordered=True)

    # Build mean-expr + frac-expr matrix per (gene, group, celltype)
    all_listed = [g for cat in GENE_CATEGORIES.values() for g in cat]
    available = [g for g in all_listed if g in sub.var_names]
    logging.info('Genes in panel: %d / %d listed', len(available), len(all_listed))

    group_col = sub.obs['group'].astype(str).values
    ct_col = sub.obs['celltype'].astype(str).values
    rows = []
    for gene in available:
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        df = pd.DataFrame({'expr': v, 'group': group_col, 'celltype': ct_col})
        agg = (df.groupby(['group', 'celltype'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda x: float((x > 0).mean())))
                 .reset_index())
        agg['gene'] = gene
        rows.append(agg)
    dot_df = pd.concat(rows, ignore_index=True)

    # ---- Filter to "interesting" genes ----
    # 1) Cell-type Δ: max across groups of |Neuron - Astrocyte|
    ct_pivot = dot_df.pivot_table(index=['gene', 'group'], columns='celltype', values='mean_expr')
    ct_pivot['delta'] = (ct_pivot['Neuron'] - ct_pivot['Astrocyte']).abs()
    ct_delta = ct_pivot.groupby('gene')['delta'].max()

    # 2) Treatment Δ: |EtOH_veh - EtOH_MCT1i| within each celltype
    treat_pivot = dot_df[dot_df['group'].isin(['EtOH_veh', 'EtOH_MCT1i'])].pivot_table(
        index=['gene', 'celltype'], columns='group', values='mean_expr')
    treat_pivot['delta'] = (treat_pivot['EtOH_veh'] - treat_pivot['EtOH_MCT1i']).abs()
    treat_delta = treat_pivot.groupby('gene')['delta'].max()

    score_df = pd.DataFrame({
        'celltype_delta': ct_delta.fillna(0.0),
        'treatment_delta': treat_delta.fillna(0.0),
    })
    score_df['pass_celltype'] = score_df['celltype_delta'] >= args.celltype_delta
    score_df['pass_treatment'] = score_df['treatment_delta'] >= args.treatment_delta
    score_df['pass_any'] = score_df['pass_celltype'] | score_df['pass_treatment']
    score_csv = args.out_dir / f'{args.label}_FocusedDotPlot_change_scores_{args.date}.csv'
    score_df.to_csv(score_csv)
    logging.info('Wrote change-score table: %s', score_csv)
    keep_genes = set(score_df[score_df['pass_any']].index)

    # Restrict to passing genes, in category order
    cat_to_genes = {cat: [g for g in gs if g in keep_genes] for cat, gs in GENE_CATEGORIES.items()}
    cat_to_genes = {k: v for k, v in cat_to_genes.items() if v}
    ordered_genes = [g for gs in cat_to_genes.values() for g in gs]
    logging.info('Kept %d / %d genes after filtering', len(ordered_genes), len(available))
    logging.info('  Pass celltype only: %d', int((score_df['pass_celltype'] & ~score_df['pass_treatment']).sum()))
    logging.info('  Pass treatment only: %d', int((~score_df['pass_celltype'] & score_df['pass_treatment']).sum()))
    logging.info('  Pass both: %d', int((score_df['pass_celltype'] & score_df['pass_treatment']).sum()))

    dot_df = dot_df[dot_df['gene'].isin(ordered_genes)].copy()
    dot_csv = args.out_dir / f'{args.label}_FocusedDotPlot_NeuronAstrocyte_{args.date}.csv'
    dot_df.to_csv(dot_csv, index=False)
    logging.info('Wrote dot plot source data: %s', dot_csv)

    # ---- Plot ----
    x_neu = list(range(len(GROUPS)))
    x_ast = [g + len(GROUPS) + 1 for g in range(len(GROUPS))]
    col_specs = [(g, 'Neuron', x_neu[i]) for i, g in enumerate(GROUPS)] + \
                [(g, 'Astrocyte', x_ast[i]) for i, g in enumerate(GROUPS)]

    n_genes = len(ordered_genes)
    fig_w = 1.6 + (2 * len(GROUPS) + 1) * 0.45 + 1.8
    fig_h = 1.5 + n_genes * 0.22
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    vmin = float(dot_df['mean_expr'].min())
    vmax = float(dot_df['mean_expr'].max())
    dot_idx = dot_df.set_index(['gene', 'group', 'celltype'])
    for i, gene in enumerate(ordered_genes):
        for grp, ct, x in col_specs:
            try:
                m = dot_idx.loc[(gene, grp, ct), 'mean_expr']
                f = dot_idx.loc[(gene, grp, ct), 'frac_expr']
            except KeyError:
                continue
            if pd.isna(m) or pd.isna(f):
                continue
            ax.scatter(x, i, s=(f * 220) + 5, c=[m], cmap='Reds',
                       vmin=vmin, vmax=vmax, edgecolor='black', linewidth=0.4)

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(ordered_genes, fontsize=8)
    ax.invert_yaxis()
    ax.set_xticks(x_neu + x_ast)
    ax.set_xticklabels(GROUPS + GROUPS, rotation=45, ha='right', fontsize=8)
    ax.set_xlim(-0.8, x_ast[-1] + 0.8)
    ax.set_ylim(n_genes - 0.5, -2.0)

    sep_x = len(GROUPS) - 0.5 + 1.0
    ax.axvline(sep_x, color='black', lw=0.8, alpha=0.5)
    ax.text((x_neu[0] + x_neu[-1]) / 2, -1.4, 'Neuron',
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text((x_ast[0] + x_ast[-1]) / 2, -1.4, 'Astrocyte',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

    label_x = x_ast[-1] + 1.0
    y_cur = 0
    cat_lines = []
    for cat, gs in cat_to_genes.items():
        n = len(gs)
        center = y_cur + (n - 1) / 2
        ax.text(label_x, center, cat, fontsize=8, va='center', ha='left', fontweight='bold')
        y_cur += n
        if y_cur < n_genes:
            cat_lines.append(y_cur - 0.5)
    for y in cat_lines:
        ax.axhline(y, color='black', lw=0.4, alpha=0.35)

    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.18, anchor=(1.0, 0.0))
    cbar.set_label('Mean log-norm expression', fontsize=8)

    size_ax = fig.add_axes([0.92, 0.55, 0.08, 0.22])
    size_ax.axis('off')
    size_ax.set_title('% cells\nexpressing', fontsize=8, loc='left')
    for k, frac in enumerate([0.1, 0.25, 0.5, 0.75, 1.0]):
        size_ax.scatter(0.2, k, s=(frac * 220) + 5, c='lightgray',
                        edgecolor='black', linewidth=0.4)
        size_ax.text(0.55, k, f'{int(frac * 100)}%', va='center', fontsize=7)
    size_ax.set_xlim(0, 1); size_ax.set_ylim(-0.5, 4.5)

    fig.suptitle(
        f'{args.label}: focused dot plot — curated genes that vary visually\n'
        f'(filtered: |N–A| ≥ {args.celltype_delta:.1f} log-norm in some group  OR  '
        f'|EtOH_veh – EtOH_MCT1i| ≥ {args.treatment_delta:.1f} log-norm in N or A)',
        fontsize=10)
    pdf_out = args.out_dir / f'{args.label}_FocusedDotPlot_NeuronAstrocyte_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', pdf_out)


if __name__ == '__main__':
    main()
