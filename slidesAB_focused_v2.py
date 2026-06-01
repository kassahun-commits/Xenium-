#!/usr/bin/env python3
"""
Focused (v2) — combined Slide A + Slide B.

Adds force-include list (genes that are always shown even if they don't
"change much") and produces TWO outputs from one run:

  1. Focused DOT plot (same look as v1) with updated gene list.
  2. Per-sample mean BAR plot grid — one panel per gene, x = 7 groups,
     side-by-side Neuron vs Astrocyte bars, error bars = SEM across mice,
     individual mouse means overlaid as black dots. Much clearer than a
     cell-level box plot when data are zero-inflated.

Both outputs use the same final gene list so the figures line up.
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
from matplotlib.patches import Patch

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
GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
          'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
CELLTYPES = ['Neuron', 'Astrocyte']
CT_COLORS = {'Neuron': '#1f77b4', 'Astrocyte': '#d62728'}

GENE_CATEGORIES = {
    'Cell-type identity': ['Rbfox3', 'Snap25', 'Syn1', 'Map2',
                           'Gfap', 'Aqp4', 'Slc1a3', 'S100b', 'Aldh1l1'],
    'Alcohol & acetaldehyde metab': ['Adh1', 'Adh5', 'Aldh2', 'Aldh1a1',
                                     'Cyp2e1', 'Cat', 'Akr1a1'],
    'Acetate / Ac-CoA metab': ['Acss1', 'Acss2', 'Acly', 'Acaca', 'Acacb'],
    'Lactate / MCT axis': ['Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb'],
    'Glucose transport': ['Slc2a1', 'Slc2a2', 'Slc2a3', 'Slc2a4', 'Slc2a5', 'Slc5a2'],
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
    'Glutamate transporters (EAAT/Slc1a)': ['Slc1a1', 'Slc1a3'],
    'HATs (lysine acetyltransferases)': ['Ep300', 'Crebbp', 'Kat2a', 'Kat2b',
                                          'Kat5', 'Kat6a', 'Kat6b', 'Kat7', 'Kat8',
                                          'Clock', 'Ncoa1', 'Ncoa2', 'Ncoa3'],
    'HDACs / Sirtuins': ['Hdac1', 'Hdac2', 'Hdac3', 'Hdac4', 'Hdac5', 'Hdac6',
                          'Hdac8', 'Hdac9',
                          'Sirt1', 'Sirt2', 'Sirt3', 'Sirt4', 'Sirt5', 'Sirt6', 'Sirt7'],
}

# Genes that are scientifically important — always include even if they
# don't change visually. Edit this list if you want to pin more.
DEFAULT_FORCE_INCLUDE = [
    'Slc16a1',  # MCT1
    'Slc16a3',  # MCT4
    'Slc16a7',  # MCT2
    'Slc2a1',   # GLUT1
    'Slc2a2',   # GLUT2
    'Slc2a3',   # GLUT3
    'Slc2a4',   # GLUT4
    'Slc2a5',   # GLUT5
    'Slc5a2',   # SGLT2
    'Acss1', 'Acss2',
    'Adh5', 'Cat',  # core alcohol enzymes
    'Mat2a', 'Mat2b', 'Ahcy',  # MAT2A axis
    # HATs - PI-relevant chromatin writers
    'Ep300', 'Crebbp', 'Kat2a', 'Kat5',
    # HDACs / Sirtuins - chromatin erasers
    'Hdac1', 'Hdac2', 'Hdac3', 'Sirt1',
    # Glutamate transporters
    'Slc1a1', 'Slc1a3',
]


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


# -------------------- Plotting helpers --------------------

def plot_focused_dotplot(dot_df, ordered_genes, cat_to_genes, out_pdf, title):
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
    ax.set_yticks(range(n_genes)); ax.set_yticklabels(ordered_genes, fontsize=8)
    ax.invert_yaxis()
    ax.set_xticks(x_neu + x_ast); ax.set_xticklabels(GROUPS + GROUPS, rotation=45, ha='right', fontsize=8)
    ax.set_xlim(-0.8, x_ast[-1] + 0.8); ax.set_ylim(n_genes - 0.5, -2.0)
    sep_x = len(GROUPS) - 0.5 + 1.0
    ax.axvline(sep_x, color='black', lw=0.8, alpha=0.5)
    ax.text((x_neu[0] + x_neu[-1]) / 2, -1.4, 'Neuron', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text((x_ast[0] + x_ast[-1]) / 2, -1.4, 'Astrocyte', ha='center', va='bottom', fontsize=12, fontweight='bold')
    label_x = x_ast[-1] + 1.0
    y_cur = 0; cat_lines = []
    for cat, gs in cat_to_genes.items():
        n = len(gs); center = y_cur + (n - 1) / 2
        ax.text(label_x, center, cat, fontsize=8, va='center', ha='left', fontweight='bold')
        y_cur += n
        if y_cur < n_genes: cat_lines.append(y_cur - 0.5)
    for y in cat_lines:
        ax.axhline(y, color='black', lw=0.4, alpha=0.35)
    for spine in ('top', 'right'): ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.18, anchor=(1.0, 0.0))
    cbar.set_label('Mean log-norm expression', fontsize=8)
    size_ax = fig.add_axes([0.92, 0.55, 0.08, 0.22])
    size_ax.axis('off')
    size_ax.set_title('% cells\nexpressing', fontsize=8, loc='left')
    for k, frac in enumerate([0.1, 0.25, 0.5, 0.75, 1.0]):
        size_ax.scatter(0.2, k, s=(frac * 220) + 5, c='lightgray', edgecolor='black', linewidth=0.4)
        size_ax.text(0.55, k, f'{int(frac * 100)}%', va='center', fontsize=7)
    size_ax.set_xlim(0, 1); size_ax.set_ylim(-0.5, 4.5)
    fig.suptitle(title, fontsize=10)
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()


def plot_persample_bargrid(sample_means_df, ordered_genes, cat_to_genes, out_pdf, title):
    """Grid of bar plots, one panel per gene.
    sample_means_df rows: gene, group, celltype, sample_id, mean_expr_log."""
    n_genes = len(ordered_genes)
    n_cols = 7
    n_rows = (n_genes + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.4 * n_cols, 1.8 * n_rows), squeeze=False)
    np.random.seed(0)

    width = 0.32
    offsets = {'Neuron': -0.20, 'Astrocyte': 0.20}

    for i, gene in enumerate(ordered_genes):
        r, c = divmod(i, n_cols)
        ax = axes[r, c]
        d_gene = sample_means_df[sample_means_df['gene'] == gene]
        if len(d_gene) == 0:
            ax.set_title(gene + ' (missing)', fontsize=8); ax.axis('off'); continue

        max_val = 0.0
        for gi, grp in enumerate(GROUPS):
            for ct in CELLTYPES:
                d = d_gene[(d_gene['group'] == grp) & (d_gene['celltype'] == ct)]
                if len(d) == 0: continue
                pos = gi + offsets[ct]
                means = d['mean_expr_log'].values
                bar_h = means.mean()
                sem = means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1 else 0
                ax.bar(pos, bar_h, width=width,
                       color=CT_COLORS[ct], edgecolor='black', linewidth=0.5,
                       alpha=0.7, zorder=2)
                if sem > 0:
                    ax.errorbar(pos, bar_h, yerr=sem, color='black',
                                capsize=2, lw=0.8, zorder=3)
                # Overlay individual mouse means
                jitter = np.random.uniform(-0.06, 0.06, size=len(means))
                ax.scatter([pos] * len(means) + jitter, means,
                           c='black', s=10, zorder=4,
                           edgecolor='white', linewidth=0.4)
                max_val = max(max_val, bar_h + sem, means.max() if len(means) else 0)

        ax.set_xticks(range(len(GROUPS)))
        ax.set_xticklabels(GROUPS, rotation=45, ha='right', fontsize=6)
        ax.set_title(gene, fontsize=9, fontweight='bold')
        ax.tick_params(axis='y', labelsize=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.6, len(GROUPS) - 0.4)
        ax.set_ylim(0, max(max_val * 1.18, 0.1))
        if c == 0:
            ax.set_ylabel('mean log-norm expr (per-mouse)', fontsize=6)

    # Hide unused
    for j in range(n_genes, n_rows * n_cols):
        r, c = divmod(j, n_cols)
        axes[r, c].axis('off')

    legend_handles = [
        Patch(facecolor=CT_COLORS['Neuron'], edgecolor='black', alpha=0.7, label='Neuron'),
        Patch(facecolor=CT_COLORS['Astrocyte'], edgecolor='black', alpha=0.7, label='Astrocyte'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                   markersize=6, label='per-mouse mean'),
        plt.Line2D([0], [0], color='black', lw=0.8, label='SEM'),
    ]
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=4, frameon=False, fontsize=10)
    fig.suptitle(title, fontsize=11, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()


# -------------------- Main --------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--celltype-delta', type=float, default=0.5)
    ap.add_argument('--treatment-delta', type=float, default=0.3)
    ap.add_argument('--force-include', nargs='*', default=DEFAULT_FORCE_INCLUDE)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.StreamHandler(sys.stdout)])

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

    # ---- Per-cell expr → per-(group,celltype) summary + per-(sample,celltype) mean ----
    all_listed = [g for cat in GENE_CATEGORIES.values() for g in cat]
    available = list(dict.fromkeys(g for g in all_listed if g in sub.var_names))  # dedup, preserve order
    logging.info('Genes in panel: %d / %d', len(available), len(all_listed))

    group_col = sub.obs['group'].astype(str).values
    ct_col = sub.obs['celltype'].astype(str).values
    sample_col = sub.obs['sample_id'].astype(str).values

    dot_rows = []
    sample_mean_rows = []
    for gene in available:
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        # group×celltype summary
        df = pd.DataFrame({'expr': v, 'group': group_col, 'celltype': ct_col, 'sample_id': sample_col})
        agg = (df.groupby(['group', 'celltype'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda x: float((x > 0).mean())))
                 .reset_index())
        agg['gene'] = gene
        dot_rows.append(agg)
        # per-sample × celltype mean
        smean = df.groupby(['group', 'celltype', 'sample_id'])['expr'].mean().reset_index()
        smean = smean.rename(columns={'expr': 'mean_expr_log'})
        smean['gene'] = gene
        sample_mean_rows.append(smean)
    dot_df = pd.concat(dot_rows, ignore_index=True)
    sample_df = pd.concat(sample_mean_rows, ignore_index=True)

    # ---- Filter ----
    ct_pivot = dot_df.pivot_table(index=['gene', 'group'], columns='celltype', values='mean_expr')
    ct_pivot['delta'] = (ct_pivot['Neuron'] - ct_pivot['Astrocyte']).abs()
    ct_delta = ct_pivot.groupby('gene')['delta'].max()
    treat_pivot = dot_df[dot_df['group'].isin(['EtOH_veh', 'EtOH_MCT1i'])].pivot_table(
        index=['gene', 'celltype'], columns='group', values='mean_expr')
    treat_pivot['delta'] = (treat_pivot['EtOH_veh'] - treat_pivot['EtOH_MCT1i']).abs()
    treat_delta = treat_pivot.groupby('gene')['delta'].max()
    score_df = pd.DataFrame({'celltype_delta': ct_delta.fillna(0.0),
                             'treatment_delta': treat_delta.fillna(0.0)})
    score_df['pass_celltype'] = score_df['celltype_delta'] >= args.celltype_delta
    score_df['pass_treatment'] = score_df['treatment_delta'] >= args.treatment_delta
    score_df['pass_any'] = score_df['pass_celltype'] | score_df['pass_treatment']
    force_set = {g for g in args.force_include if g in available}
    keep_genes = set(score_df[score_df['pass_any']].index) | force_set
    logging.info('Force-included (present in panel): %s', sorted(force_set))
    logging.info('Force-included but NOT in panel: %s',
                 sorted(set(args.force_include) - force_set))

    # Dedup: genes can appear in multiple categories; keep first occurrence
    seen: set = set()
    cat_to_genes = {}
    for cat, gs in GENE_CATEGORIES.items():
        keep = [g for g in gs if g in keep_genes and g not in seen]
        seen.update(keep)
        if keep:
            cat_to_genes[cat] = keep
    ordered_genes = [g for gs in cat_to_genes.values() for g in gs]
    logging.info('Total genes kept: %d', len(ordered_genes))

    # Save source data
    dot_df_kept = dot_df[dot_df['gene'].isin(ordered_genes)].copy()
    sample_df_kept = sample_df[sample_df['gene'].isin(ordered_genes)].copy()
    dot_df_kept.to_csv(args.out_dir / f'{args.label}_FocusedV2_dotplot_data_{args.date}.csv', index=False)
    sample_df_kept.to_csv(args.out_dir / f'{args.label}_FocusedV2_persample_means_{args.date}.csv', index=False)
    score_df.to_csv(args.out_dir / f'{args.label}_FocusedV2_change_scores_{args.date}.csv')

    # ---- Plot dot plot ----
    title_dot = (f'{args.label}: focused dot plot v2 — varying genes + key alcohol/metab axis\n'
                 f'(filter ΔN-A≥{args.celltype_delta:.1f} OR ΔEtOH_veh-EtOH_MCT1i≥{args.treatment_delta:.1f}; '
                 f'plus {len(force_set)} force-included)')
    plot_focused_dotplot(
        dot_df_kept, ordered_genes, cat_to_genes,
        out_pdf=args.out_dir / f'{args.label}_FocusedDotPlot_v2_NeuronAstrocyte_{args.date}.pdf',
        title=title_dot)

    # ---- Plot per-sample bar grid ----
    title_bar = (f'{args.label}: per-mouse mean bar plots — same gene list as dot plot v2\n'
                 f'(bars = mean of per-mouse means; error bars = SEM; dots = individual mice; n=2-3 per group)')
    plot_persample_bargrid(
        sample_df_kept, ordered_genes, cat_to_genes,
        out_pdf=args.out_dir / f'{args.label}_FocusedBarPlot_v2_NeuronAstrocyte_{args.date}.pdf',
        title=title_bar)

    logging.info('Done.')


if __name__ == '__main__':
    main()
