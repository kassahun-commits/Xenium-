#!/usr/bin/env python3
"""
Combined Slide A + Slide B curated dot plot.

Loads both slides' counts + Xenium-Explorer ROI exports, assigns cells
to ROIs (point-in-polygon), maps ROI labels -> (group, replicate),
concatenates into a single AnnData (with slide annotation), runs the
marker-based cell typing on the combined data, and produces the curated
~100-gene dot plot in Neuron vs Astrocyte blocks across 7 groups, now
with proper biological replicates (n=2-3 per group).

APP is not included (no APP ROIs drawn on Slide A by user choice).
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

# ----------------------------- Config -----------------------------

LABEL_MAP = {
    # --- Slide B labels ---
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
    # --- Slide A labels ---
    'EtOH +MCT1i 2':  ('EtOH_MCT1i', 2),
    'EtOH +MCT1i 3':  ('EtOH_MCT1i', 3),
    'H20+veh 2':      ('H2O_veh',    2),
    'H20 +MCT1i 2':   ('H2O_MCT1i',  2),
    'Chronic EtOH 2': ('ChronicEtOH', 2),
    'chronic EtOH 3': ('ChronicEtOH', 3),  # lowercase 'c' in user label
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
    'Cell-type identity': [
        'Rbfox3', 'Snap25', 'Syn1', 'Map2',
        'Gfap', 'Aqp4', 'Slc1a3', 'S100b', 'Aldh1l1',
    ],
    'Alcohol & acetaldehyde metab': [
        'Adh1', 'Adh5', 'Aldh2', 'Aldh1a1', 'Cyp2e1', 'Cat', 'Akr1a1',
    ],
    'Acetate / Ac-CoA metab': [
        'Acss1', 'Acss2', 'Acly', 'Acaca', 'Acacb',
    ],
    'Lactate / MCT axis': [
        'Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb',
    ],
    'Glucose transport': [
        'Slc2a1', 'Slc2a2', 'Slc2a3', 'Slc2a4', 'Slc2a5', 'Slc5a2',
    ],
    'Glycolysis / TCA': [
        'Gapdh', 'Pgk1', 'Pkm', 'Hk1', 'Hk2', 'Pfkm', 'Aldoa',
        'Eno1', 'Eno2', 'Idh1', 'Idh2', 'Cs',
    ],
    'SAM / MAT2A axis': [
        'Mat1a', 'Mat2a', 'Mat2b', 'Ahcy', 'Mthfr', 'Cbs',
    ],
    'Glutathione / oxidative stress': [
        'Gclc', 'Gclm', 'Gpx1', 'Gpx4', 'Gss', 'Gstm1',
        'Sod1', 'Sod2', 'Hmox1', 'Nfe2l2',
    ],
    'Stress / UPR': [
        'Atf3', 'Atf4', 'Ddit3', 'Hspa1a', 'Hsp90aa1', 'Hspb1', 'Trp53',
    ],
    'Activity / IEGs': [
        'Fos', 'Jun', 'Junb', 'Arc', 'Egr1', 'Egr2',
        'Npas4', 'Homer1', 'Bdnf',
    ],
    'Glutamate signaling': [
        'Slc17a7', 'Slc17a6', 'Slc1a2',
        'Grin1', 'Grin2a', 'Grin2b',
        'Gria1', 'Gria2',
        'Grm1', 'Grm2', 'Grm5',
        'Camk2a', 'Camk2b',
    ],
    'GABA signaling': [
        'Gad1', 'Gad2', 'Gabra1', 'Gabra2', 'Gabbr1', 'Slc32a1',
    ],
    'Endocannabinoid': [
        'Cnr1', 'Faah', 'Mgll', 'Daglb',
    ],
    'Neuroinflammation': [
        'Tnf', 'Il1b', 'Il6', 'Il10', 'Tlr4', 'Nfkb1', 'Stat3',
    ],
    'Glutamate transporters (EAAT/Slc1a)': ['Slc1a1', 'Slc1a3'],
    'HATs (lysine acetyltransferases)': [
        'Ep300', 'Crebbp', 'Kat2a', 'Kat2b', 'Kat5', 'Kat6a', 'Kat6b',
        'Kat7', 'Kat8', 'Clock', 'Ncoa1', 'Ncoa2', 'Ncoa3',
    ],
    'HDACs / Sirtuins': [
        'Hdac1', 'Hdac2', 'Hdac3', 'Hdac4', 'Hdac5', 'Hdac6', 'Hdac8', 'Hdac9',
        'Sirt1', 'Sirt2', 'Sirt3', 'Sirt4', 'Sirt5', 'Sirt6', 'Sirt7',
    ],
}


# ----------------------------- Helpers -----------------------------

def load_polygons(csv_path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(csv_path, comment='#')
    if not {'Selection', 'X', 'Y'}.issubset(df.columns):
        raise ValueError(f'Annotations CSV missing cols; got {list(df.columns)}')
    polygons = {}
    for name, sub in df.groupby('Selection', sort=False):
        polygons[name] = sub[['X', 'Y']].values.astype(float)
    return polygons


def assign_cells_to_rois(cells_df, polygons):
    out = pd.Series(np.full(len(cells_df), np.nan, dtype=object),
                    index=cells_df.index, name='roi')
    pts = cells_df[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        path = MplPath(xy)
        inside = path.contains_points(pts)
        mask = inside & ~out.notna().values
        out.values[mask] = name
    return out


def process_slide(slide_dir: Path, ann_csv: Path, slide_label: str):
    """Return AnnData for one slide with obs containing roi/group/replicate/sample_id/slide."""
    cells = pd.read_csv(slide_dir / 'cells.csv.gz').set_index('cell_id')
    polygons = load_polygons(ann_csv)
    cells['roi'] = assign_cells_to_rois(cells, polygons)
    n_in = cells['roi'].notna().sum()
    logging.info('%s: %s cells loaded; %s assigned to a ROI (%.1f%%)',
                 slide_label, f'{len(cells):,}', f'{n_in:,}',
                 100.0 * n_in / len(cells))
    # Warn about unmapped ROIs
    unmapped = set(cells['roi'].dropna().unique()) - set(LABEL_MAP.keys())
    if unmapped:
        logging.warning('%s: unmapped ROI labels (cells excluded): %s', slide_label, unmapped)
    cells['group'] = cells['roi'].map(lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['slide'] = slide_label
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)

    cells_in = cells.dropna(subset=['group']).copy()
    logging.info('%s: %s cells with a known group', slide_label, f'{len(cells_in):,}')

    # Load counts
    adata = sc.read_10x_h5(str(slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    # Prefix cell IDs with slide label to avoid collisions on concat
    adata.obs_names = [f'{slide_label}::{x}' for x in adata.obs_names]
    cells_in.index = [f'{slide_label}::{x}' for x in cells_in.index]

    common = adata.obs_names.intersection(cells_in.index)
    adata = adata[common].copy()
    adata.obs = cells_in.loc[common, ['roi', 'group', 'replicate',
                                       'sample_id', 'slide',
                                       'x_centroid', 'y_centroid']].copy()
    return adata


# ----------------------------- Main -----------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--slide-a-dir', required=True, type=Path)
    p.add_argument('--slide-a-ann', required=True, type=Path)
    p.add_argument('--slide-b-dir', required=True, type=Path)
    p.add_argument('--slide-b-ann', required=True, type=Path)
    p.add_argument('--out-dir', required=True, type=Path)
    p.add_argument('--label', default='SlidesAB')
    p.add_argument('--min-counts', type=int, default=10)
    p.add_argument('--min-cells', type=int, default=5)
    p.add_argument('--date', default=date.today().isoformat())
    return p.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== %s combined analysis ===', args.label)

    # Process each slide
    ada_a = process_slide(args.slide_a_dir, args.slide_a_ann, 'SlideA')
    ada_b = process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    logging.info('SlideA: %s cells x %s genes', f'{ada_a.shape[0]:,}', f'{ada_a.shape[1]:,}')
    logging.info('SlideB: %s cells x %s genes', f'{ada_b.shape[0]:,}', f'{ada_b.shape[1]:,}')

    # Concat (same panel -> inner join on var = same gene set)
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same',
                      label='slide_src', keys=['SlideA', 'SlideB'])
    logging.info('Combined: %s cells x %s genes', f'{adata.shape[0]:,}', f'{adata.shape[1]:,}')

    # Sample composition
    logging.info('Samples per group (after combining):')
    samp_per_group = adata.obs.groupby('group', observed=True)['sample_id'].nunique()
    logging.info('\n%s', samp_per_group.to_string())

    # QC
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    n0 = adata.shape[0]; g0 = adata.shape[1]
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    logging.info('After QC: %s cells (dropped %s) x %s genes (dropped %s)',
                 f'{adata.shape[0]:,}', f'{n0 - adata.shape[0]:,}',
                 f'{adata.shape[1]:,}', f'{g0 - adata.shape[1]:,}')

    # Normalize + log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Cell-type scoring
    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gs_present,
                          score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    score_cols = [f'score_{ct}' for ct in MARKERS]
    scores = adata.obs[score_cols]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_val = scores.max(axis=1)
    best_ct[best_val <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(
        best_ct, categories=list(MARKERS.keys()) + ['Unclassified'])

    # Counts per sample × celltype (sanity)
    ct_counts = (adata.obs.groupby(['slide', 'group', 'replicate', 'celltype'],
                                   observed=True).size().unstack(fill_value=0))
    logging.info('Cell type counts per sample:\n%s', ct_counts.to_string())
    ct_counts.to_csv(args.out_dir / f'{args.label}_celltype_counts_per_sample_{args.date}.csv')

    # Per-cell CSV
    obs_out = adata.obs.copy()
    obs_out.to_csv(args.out_dir / f'{args.label}_cells_with_roi_{args.date}.csv')

    # Restrict to Neuron + Astrocyte for dot plot
    sub = adata[adata.obs['celltype'].astype(str).isin(CELLTYPES)].copy()
    sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str), categories=GROUPS, ordered=True)
    logging.info('Cells for dot plot (N+A in 7 groups): %s', f'{sub.shape[0]:,}')

    # Curated gene list resolution
    all_listed = [g for cat in GENE_CATEGORIES.values() for g in cat]
    available = list(dict.fromkeys(g for g in all_listed if g in sub.var_names))  # dedup, preserve order
    missing = [g for g in all_listed if g not in sub.var_names]
    logging.info('Genes listed: %d; available: %d; missing %d: %s',
                 len(all_listed), len(available), len(missing), missing)

    # Dedup: genes can appear in multiple categories; keep first occurrence
    seen: set = set()
    cat_to_genes = {}
    for cat, gs in GENE_CATEGORIES.items():
        keep = [g for g in gs if g in available and g not in seen]
        seen.update(keep)
        if keep:
            cat_to_genes[cat] = keep
    ordered_genes = [g for gs in cat_to_genes.values() for g in gs]

    # Aggregate per (group, celltype, gene)
    group_col = sub.obs['group'].astype(str).values
    ct_col = sub.obs['celltype'].astype(str).values
    rows = []
    for gene in ordered_genes:
        X = sub[:, gene].X
        v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
        df = pd.DataFrame({'expr': v, 'group': group_col, 'celltype': ct_col})
        agg = (df.groupby(['group', 'celltype'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda x: float((x > 0).mean())),
                      n_cells=('expr', 'size'),
                      n_expr=('expr', lambda x: int((x > 0).sum())))
                 .reset_index())
        agg['gene'] = gene
        rows.append(agg)
    dot_df = pd.concat(rows, ignore_index=True)
    dot_csv = args.out_dir / f'{args.label}_CuratedDotPlot_NeuronAstrocyte_{args.date}.csv'
    dot_df.to_csv(dot_csv, index=False)
    logging.info('Wrote %s', dot_csv)

    # ----- Plot -----
    x_neu = list(range(len(GROUPS)))
    x_ast = [g + len(GROUPS) + 1 for g in range(len(GROUPS))]
    col_specs = [(g, 'Neuron', x_neu[i]) for i, g in enumerate(GROUPS)] + \
                [(g, 'Astrocyte', x_ast[i]) for i, g in enumerate(GROUPS)]

    n_genes = len(ordered_genes)
    fig_w = 1.6 + (2 * len(GROUPS) + 1) * 0.45 + 1.8
    fig_h = 1.5 + n_genes * 0.20
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
                       vmin=vmin, vmax=vmax,
                       edgecolor='black', linewidth=0.4)

    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(ordered_genes, fontsize=7)
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
        ax.text(label_x, center, cat, fontsize=8, va='center',
                ha='left', fontweight='bold')
        y_cur += n
        if y_cur < n_genes:
            cat_lines.append(y_cur - 0.5)
    for y in cat_lines:
        ax.axhline(y, color='black', lw=0.4, alpha=0.35)

    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis='y', length=0)

    sm = plt.cm.ScalarMappable(cmap='Reds',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.18, anchor=(1.0, 0.0))
    cbar.set_label('Mean log-norm expression', fontsize=8)

    size_ax = fig.add_axes([0.92, 0.65, 0.08, 0.18])
    size_ax.axis('off')
    size_ax.set_title('% cells\nexpressing', fontsize=8, loc='left')
    for k, frac in enumerate([0.1, 0.25, 0.5, 0.75, 1.0]):
        size_ax.scatter(0.2, k, s=(frac * 220) + 5, c='lightgray',
                        edgecolor='black', linewidth=0.4)
        size_ax.text(0.55, k, f'{int(frac * 100)}%', va='center', fontsize=7)
    size_ax.set_xlim(0, 1); size_ax.set_ylim(-0.5, 4.5)

    # Build replicate-count annotation under group labels
    rep_counts = (sub.obs[['group', 'sample_id']].drop_duplicates()
                  .groupby('group', observed=True).size().to_dict())
    rep_str = lambda g: f'n={rep_counts.get(g, 0)}'
    for x, g in zip(x_neu, GROUPS):
        ax.text(x, n_genes + 0.2, rep_str(g), ha='center', va='top',
                fontsize=6, color='#444')
    for x, g in zip(x_ast, GROUPS):
        ax.text(x, n_genes + 0.2, rep_str(g), ha='center', va='top',
                fontsize=6, color='#444')

    fig.suptitle(
        f'{args.label}: curated alcohol / metabolism / activity gene panel '
        f'across conditions in Neurons vs Astrocytes\n'
        f'(combined Slide A + Slide B; n shown below each column)',
        fontsize=11)
    pdf_out = args.out_dir / f'{args.label}_CuratedDotPlot_NeuronAstrocyte_{args.date}.pdf'
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()
    logging.info('Wrote %s', pdf_out)
    logging.info('Done.')


if __name__ == '__main__':
    main()
