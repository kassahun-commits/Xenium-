#!/usr/bin/env python3
"""
Multi-page per-mouse mean BAR PLOTS for excitatory & inhibitory neurons across
all 7 condition groups, using a custom gene list.

Differences vs prior box-plot version:
  * BAR plots (height = mean of per-mouse means; error bar = SEM; black dots
    = individual mouse means) — cleaner for zero-inflated spatial data.
  * All 7 condition groups on the x-axis (not just acute/chronic):
    H2O_veh, H2O_MCT1i, EtOH_veh, EtOH_MCT1i, ChronicEtOH, MAT2A_CM, MAT2A_OE.
  * Multi-page PDF: 20 genes per page (5 cols x 4 rows).
  * Custom gene list provided by user — plus cell-type / synaptic markers
    and GFP (MAT2A_OE reporter).
  * Gene labels use 'gene (alias)' format where applicable.

Outputs (per subtype):
  - SlidesAB_<Subtype>_AllGroups_BarPlots_<date>.pdf  (multi-page)
  - SlidesAB_<Subtype>_AllGroups_persample_means_<date>.csv
  - SlidesAB_<Subtype>_log.txt
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
from matplotlib.patches import Patch
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
INHIB_MARKERS = ['Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln', 'Lhx6']

ALL_GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
              'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE']
GROUP_COLORS = {
    'H2O_veh':     '#888888',
    'H2O_MCT1i':   '#17becf',
    'EtOH_veh':    '#ff7f0e',
    'EtOH_MCT1i':  '#2ca02c',
    'ChronicEtOH': '#a83232',
    'MAT2A_CM':    '#9467bd',
    'MAT2A_OE':    '#8c564b',
}

ALIAS = {
    'Slc16a1': 'MCT1', 'Slc16a7': 'MCT2', 'Slc16a3': 'MCT4',
    'Slc17a7': 'VGLUT1', 'Slc17a6': 'VGLUT2',
    'Slc2a1': 'GLUT1', 'Slc2a3': 'GLUT3', 'Slc2a4': 'GLUT4',
    'Slc1a3': 'GLAST', 'Slc1a1': 'EAAC1',
    'Slc32a1': 'VGAT', 'Rbfox3': 'NeuN',
    'Crebbp': 'CBP', 'Ep300': 'p300',
    'Kat2a': 'GCN5', 'Kat2b': 'PCAF', 'Kat5': 'Tip60',
    'GFP': 'reporter',
}

# User-provided gene list (in given order)
USER_GENES = [
    'Ahcy', 'Atf3', 'Brd9', 'Btg2', 'Chd3', 'Chd8', 'Chdh', 'Dpf1', 'Fos', 'Fosb',
    'Gadd45b', 'Gnmt', 'Ino80', 'Jun', 'Junb', 'Kat7', 'Kmt5a', 'Mat1a', 'Mat2a',
    'Mat2b', 'Mbd1', 'Mbd4', 'Mtap', 'Mthfd1', 'Mthfd2', 'Mtr', 'Mtrr', 'Nsd2',
    'Nsd3', 'Pemt', 'Phgdh', 'Prdm2', 'Prdm8', 'Prmt2', 'Prmt6', 'Prmt8', 'Prmt9',
    'Psph', 'Setdb2', 'Shmt2', 'Slc25a32', 'Ss18', 'Suv39h2', 'Tyms',
    'Grik5', 'Rgs14', 'Pcp4', 'Cacng5', 'Cacnb4', 'Ccnt2', 'Cdk2ap1', 'Celf2',
    'Dlg2', 'Hmgb3', 'Hmgn3', 'Ilf3', 'Kcnh7', 'Lrrtm2', 'Lrrtm4', 'Malat1',
    'Mbnl2', 'Pclo', 'Pnisr', 'Prpf39', 'Rbm34', 'Rsrc2', 'Slc25a40', 'Tia1',
    'Zranb2',
    'Slc5a8', 'Emb', 'Bsg', 'Slc4a4', 'Crat', 'Crot', 'Slc10a2', 'Acat1', 'Acat2',
    'Slc16a7', 'Gapdh', 'Slc13a5', 'Pdha1', 'Pdhb', 'Gpx4', 'Gpx1', 'Glul', 'Adh1',
    'Unc13a', 'Abcc3', 'Abhd14b', 'Tagln3', 'Nptx1', 'Pfkm', 'Cfhr1', 'Dgkb',
    'Gng12', 'Elovl5',
]
# Cell-type identity markers prepended for context (helpful for sanity-checking
# which cells we're looking at). Then GFP, then user list.
GENE_LIST = (
    ['Rbfox3', 'Snap25', 'Syn1', 'Map2',
     'Slc17a7', 'Slc17a6', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6',
     'Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln',
     'Gfap', 'Aqp4', 'Slc1a3', 'Slc1a1',
     'Slc16a1', 'Slc16a7', 'Slc16a3',
     'GFP',]
    + [g for g in USER_GENES if g not in
       {'Rbfox3', 'Snap25', 'Syn1', 'Map2',
        'Slc17a7', 'Slc17a6', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6',
        'Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln',
        'Gfap', 'Aqp4', 'Slc1a3', 'Slc1a1',
        'Slc16a1', 'Slc16a7', 'Slc16a3',
        'GFP'}]
)


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


def plot_barpage(genes_page, sample_df, page_idx, n_pages, subtype, out_pdf_page):
    """Render one page of bar plots into the open PdfPages."""
    n_genes = len(genes_page)
    n_cols = 5; n_rows = 4
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 2.8, n_rows * 2.2),
                              squeeze=False)
    np.random.seed(0 + page_idx)
    for i, gene in enumerate(genes_page):
        r_, c_ = divmod(i, n_cols)
        ax = axes[r_, c_]
        d = sample_df[sample_df['gene'] == gene]
        if len(d) == 0:
            ax.set_title(disp(gene) + ' (missing)', fontsize=8); ax.axis('off'); continue

        offsets = {g: idx for idx, g in enumerate(ALL_GROUPS)}
        y_max = 0.0
        for grp in ALL_GROUPS:
            dg = d[d['group'] == grp]
            if len(dg) == 0: continue
            pos = offsets[grp]
            means = dg['mean_expr_log'].values
            mh = means.mean()
            sem = means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1 else 0.0
            ax.bar(pos, mh, width=0.7, color=GROUP_COLORS[grp],
                   edgecolor='black', linewidth=0.5, alpha=0.78, zorder=2)
            if sem > 0:
                ax.errorbar(pos, mh, yerr=sem, color='black', capsize=2, lw=0.8, zorder=3)
            jitter = np.random.uniform(-0.10, 0.10, size=len(means))
            ax.scatter([pos] * len(means) + jitter, means,
                       c='black', s=9, zorder=4, edgecolor='white', linewidth=0.4)
            y_max = max(y_max, mh + sem, means.max() if len(means) else 0)

        ax.set_xticks(range(len(ALL_GROUPS)))
        ax.set_xticklabels(ALL_GROUPS, rotation=45, ha='right', fontsize=6)
        ax.set_title(disp(gene), fontsize=9, fontweight='bold')
        ax.tick_params(axis='y', labelsize=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.6, len(ALL_GROUPS) - 0.4)
        ax.set_ylim(0, max(y_max * 1.20, 0.05))
        if c_ == 0:
            ax.set_ylabel('mean log-norm expr\n(per-mouse)', fontsize=6)

    for j in range(n_genes, n_rows * n_cols):
        r_, c_ = divmod(j, n_cols)
        axes[r_, c_].axis('off')

    legend_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor='black',
                            alpha=0.78, label=g) for g in ALL_GROUPS]
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                      markerfacecolor='black', markersize=6,
                                      label='per-mouse mean'))
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=0.8, label='SEM'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=5, frameon=False, fontsize=8)
    fig.suptitle(f'{subtype} neurons — per-mouse mean bar plots — '
                 f'page {page_idx + 1} of {n_pages}',
                 fontsize=11, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    out_pdf_page.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--per-page', type=int, default=20)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_SubtypeBarPlots_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Multi-page bar plots: Excitatory + Inhibitory neurons ===')

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

    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory',
                      random_state=0, n_bins=25)
    inhib_present = [g for g in INHIB_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=inhib_present, score_name='score_Inhibitory',
                      random_state=0, n_bins=25)
    is_excit = (adata.obs['celltype'].astype(str) == 'Neuron') & (adata.obs['score_Excitatory'] > 0)
    is_inhib = (adata.obs['celltype'].astype(str) == 'Neuron') & (adata.obs['score_Inhibitory'] > 0)

    # Resolve gene list against panel
    available = [g for g in GENE_LIST if g in adata.var_names]
    missing = [g for g in GENE_LIST if g not in adata.var_names]
    logging.info('Genes to plot: %d available / %d listed', len(available), len(GENE_LIST))
    if missing:
        logging.warning('Missing from panel: %s', missing)

    for subtype, mask in [('Excitatory', is_excit), ('Inhibitory', is_inhib)]:
        sub = adata[mask].copy()
        sub = sub[sub.obs['group'].astype(str).isin(ALL_GROUPS)].copy()
        sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str),
                                          categories=ALL_GROUPS, ordered=True)
        logging.info('%s neurons (after filtering to 7 groups): %d cells', subtype, sub.shape[0])
        # Counts per sample
        s_counts = sub.obs.groupby(['group', 'sample_id'], observed=True).size().reset_index(name='n')
        n_mice = {g: int(sub.obs.loc[sub.obs['group'] == g, 'sample_id'].nunique()) for g in ALL_GROUPS}
        logging.info('%s mice per group: %s', subtype, n_mice)

        # Per-sample means per (gene, group, sample_id)
        rows = []
        g_arr = sub.obs['group'].astype(str).values
        s_arr = sub.obs['sample_id'].astype(str).values
        for gene in available:
            X = sub[:, gene].X
            v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
            df = pd.DataFrame({'expr': v, 'group': g_arr, 'sample_id': s_arr})
            mean = df.groupby(['group', 'sample_id'])['expr'].mean().reset_index()
            mean.rename(columns={'expr': 'mean_expr_log'}, inplace=True)
            mean['gene'] = gene
            rows.append(mean)
        sample_df = pd.concat(rows, ignore_index=True)
        sample_csv = args.out_dir / f'{args.label}_{subtype}_AllGroups_persample_means_{args.date}.csv'
        sample_df.to_csv(sample_csv, index=False)
        logging.info('Wrote %s', sample_csv)

        # Multi-page PDF
        per_page = args.per_page
        n_pages = (len(available) + per_page - 1) // per_page
        pdf_out = args.out_dir / f'{args.label}_{subtype}_AllGroups_BarPlots_{args.date}.pdf'
        with PdfPages(pdf_out) as pdf:
            for page_idx in range(n_pages):
                chunk = available[page_idx * per_page:(page_idx + 1) * per_page]
                plot_barpage(chunk, sample_df, page_idx, n_pages, subtype, pdf)
        logging.info('Wrote %s (%d pages)', pdf_out, n_pages)


if __name__ == '__main__':
    main()
