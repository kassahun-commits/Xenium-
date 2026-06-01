#!/usr/bin/env python3
"""
Three-cell-type side-by-side BAR PLOTS for direct cross-cell-type comparison.

For each gene: one row with three small panels — Excitatory | Inhibitory |
Astrocyte — each showing per-mouse mean ± SEM across all 7 condition groups.
Shared y-axis within each row so cell-type comparisons are direct.

6 genes per page (6 rows × 3 panels = 18 panels per page).

Gene set:
  * Cell-type identity markers (for sanity)
  * MCT family + GLUT family + GFP
  * Alcohol-metabolizing enzymes
  * Acetate / Ac-CoA / glycolysis / TCA / lactate
  * Oxidative stress / detox
  * NEW: hormone receptors (GR, MR, ER, AR, thyroid, insulin/IGF, leptin, CRH, OXT, AVP)
  * NEW: metabolic/ER stress / UPR / ISR (Atf4, Atf6, Ddit3, Hspa5, Eif2ak2/3/4,
         Eif2s1, Sesn1/2, Tsc1/2, mTOR, Rheb, Foxo1/3/4, Pmaip1, Bbc3, Trib3)
  * NEW: additional oxidative stress (Keap1, Nfe2l1, Sirt3/5, Gstp1)
  * All 95 user-listed custom genes
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
SUBTYPE_LABELS = ['Excitatory', 'Inhibitory', 'Astrocyte']

ALIAS = {
    'Slc16a1': 'MCT1', 'Slc16a7': 'MCT2', 'Slc16a3': 'MCT4',
    'Slc17a7': 'VGLUT1', 'Slc17a6': 'VGLUT2',
    'Slc2a1': 'GLUT1', 'Slc2a2': 'GLUT2', 'Slc2a3': 'GLUT3',
    'Slc2a4': 'GLUT4', 'Slc2a5': 'GLUT5', 'Slc5a2': 'SGLT2',
    'Slc1a3': 'GLAST', 'Slc1a1': 'EAAC1', 'Slc1a2': 'GLT-1',
    'Slc32a1': 'VGAT', 'Slc6a11': 'GAT-3', 'Rbfox3': 'NeuN',
    'Crebbp': 'CBP', 'Ep300': 'p300',
    'Kat2a': 'GCN5', 'Kat2b': 'PCAF', 'Kat5': 'Tip60',
    'GFP': 'reporter',
    'Cat': 'catalase', 'Glul': 'glutamine synthetase',
    # Hormone receptors
    'Nr3c1': 'GR', 'Nr3c2': 'MR', 'Esr1': 'ERα', 'Esr2': 'ERβ',
    'Ar': 'AR', 'Pgr': 'PR', 'Thra': 'TRα', 'Thrb': 'TRβ',
    'Insr': 'INSR', 'Igf1r': 'IGF1R', 'Igf2r': 'IGF2R',
    'Lepr': 'LEPR', 'Ghr': 'GHR',
    'Crhr1': 'CRHR1', 'Crhr2': 'CRHR2', 'Oxtr': 'OXTR',
    'Avpr1b': 'V1bR', 'Mc4r': 'MC4R',
    # Stress
    'Atf4': 'ATF4', 'Atf6': 'ATF6', 'Ddit3': 'CHOP', 'Hspa5': 'BiP',
    'Eif2ak3': 'PERK', 'Eif2ak2': 'PKR', 'Eif2ak4': 'GCN2',
    'Eif2s1': 'eIF2α',
    'Tsc1': 'TSC1', 'Tsc2': 'TSC2', 'Mtor': 'mTOR',
    'Pmaip1': 'NOXA', 'Bbc3': 'PUMA',
}


def disp(gene: str) -> str:
    return f'{gene} ({ALIAS[gene]})' if gene in ALIAS else gene


# ---- Comprehensive gene list ----

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

# Ordered gene list, organized by theme for readability
GENE_LIST = (
    # 1. Neuronal identity / contrast
    ['Rbfox3', 'Snap25', 'Syn1', 'Map2',
     'Slc17a7', 'Slc17a6', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6',
     'Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln',
     # 2. Astrocyte identity
     'Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1',
     'Slc1a1', 'Slc6a11',
     # 3. MCT family
     'Slc16a1', 'Slc16a7', 'Slc16a3',
     # 4. Glucose transporters (GLUT1, GLUT3 emphasized per request)
     'Slc2a1', 'Slc2a3', 'Slc2a4', 'Slc2a2', 'Slc2a5', 'Slc5a2',
     # 5. Reporter
     'GFP',
     # 6. Alcohol-metabolizing enzymes
     'Adh1', 'Adh5', 'Aldh2', 'Aldh1a2', 'Cyp2e1', 'Cat',
     # 7. Acetate / Ac-CoA / glycolysis / TCA / lactate
     'Acss1', 'Acss2', 'Acly', 'Eno1', 'Eno2', 'Ldha', 'Ldhb',
     # 8. Oxidative stress / detox
     'Sod1', 'Sod2', 'Gpx1', 'Gpx4', 'Nfe2l2', 'Nfe2l1', 'Hmox1',
     'Keap1', 'Sirt3', 'Sirt5', 'Gstp1',
     # 9. Hormone receptors (NEW)
     'Nr3c1', 'Nr3c2',
     'Esr1', 'Esr2', 'Ar', 'Pgr', 'Thra', 'Thrb',
     'Insr', 'Igf1r', 'Igf2r', 'Lepr', 'Ghr',
     'Crhr1', 'Crhr2', 'Oxtr', 'Avpr1b', 'Mc4r',
     # 10. Metabolic / ER stress / UPR / ISR (NEW)
     'Atf4', 'Atf6', 'Ddit3', 'Hspa5',
     'Eif2ak3', 'Eif2ak2', 'Eif2ak4', 'Eif2s1',
     'Sesn1', 'Sesn2',
     'Tsc1', 'Tsc2', 'Mtor', 'Rheb',
     'Foxo1', 'Foxo3', 'Foxo4',
     'Pmaip1', 'Bbc3', 'Trib3',
     # 11. NEW Immediate-early genes (IEGs) — comprehensive
     # Fos family
     'Fos', 'Fosb', 'Fosl1', 'Fosl2',
     # Jun family
     'Jun', 'Junb',
     # Egr family
     'Egr1', 'Egr2', 'Egr3', 'Egr4',
     # Nr4a family (orphan nuclear receptors)
     'Nr4a1', 'Nr4a2', 'Nr4a3',
     # Activity-regulated IEGs
     'Arc', 'Bdnf', 'Homer1', 'Npas4',
     'Per1', 'Per2',
     'Dusp1', 'Dusp6',
     'Plk2', 'Plk3',
     'Nptx1', 'Nptx2',
     'Trib2', 'Trib3',
     # Other IEG-related
     'Atf3', 'Btg2', 'Gadd45b', 'Bhlhe40', 'Klf10', 'Ccn1', 'Rasgrf1',
     # 12. NEW Neurotrophic factors + receptors
     'Ntf3', 'Ntf5', 'Ngf', 'Gdnf', 'Cntf', 'Lif', 'Mst1', 'Mstn',
     'Igf1', 'Igf2',
     'Ntrk1', 'Ntrk2', 'Ntrk3', 'Ngfr',
     'Gfra1', 'Gfra2', 'Gfra3', 'Ret', 'Lifr',
    ]
    # 11. All 95 user-listed custom genes (in the order they were given)
    + [g for g in USER_GENES]
)
# Dedup while preserving order
seen = set(); GENE_LIST = [g for g in GENE_LIST if not (g in seen or seen.add(g))]


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


def render_gene_row(axes_triplet, gene, persample_by_subtype, shared_ymax=None):
    """Render 3 panels for one gene across the 3 cell types.
    persample_by_subtype: dict subtype -> DataFrame with (gene, group, sample_id, mean_expr_log).
    axes_triplet: list of 3 axes objects.
    Returns the new shared ymax.
    """
    y_max = 0.0
    for ax, subtype in zip(axes_triplet, SUBTYPE_LABELS):
        d = persample_by_subtype[subtype]
        d = d[d['gene'] == gene]
        if len(d) == 0:
            ax.text(0.5, 0.5, '(no data)', ha='center', va='center',
                    transform=ax.transAxes, fontsize=7, color='gray')
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ('top','right','left','bottom'): ax.spines[sp].set_visible(False)
            continue
        for i, grp in enumerate(ALL_GROUPS):
            dg = d[d['group'] == grp]
            if len(dg) == 0: continue
            means = dg['mean_expr_log'].values
            mh = means.mean()
            sem = means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1 else 0.0
            ax.bar(i, mh, width=0.7, color=GROUP_COLORS[grp],
                   edgecolor='black', linewidth=0.4, alpha=0.8, zorder=2)
            if sem > 0:
                ax.errorbar(i, mh, yerr=sem, color='black', capsize=1.5, lw=0.6, zorder=3)
            np.random.seed(hash((gene, subtype, grp)) % (2**32))
            jitter = np.random.uniform(-0.10, 0.10, size=len(means))
            ax.scatter([i] * len(means) + jitter, means,
                       c='black', s=6, zorder=4, edgecolor='white', linewidth=0.3)
            y_max = max(y_max, mh + sem, means.max() if len(means) else 0)
        ax.set_xticks(range(len(ALL_GROUPS)))
        ax.set_xticklabels(ALL_GROUPS, rotation=45, ha='right', fontsize=5.5)
        ax.tick_params(axis='y', labelsize=6)
        for sp in ('top','right'): ax.spines[sp].set_visible(False)
        ax.set_xlim(-0.6, len(ALL_GROUPS) - 0.4)

    # Apply shared y limits within this row
    use_max = shared_ymax if shared_ymax is not None else y_max
    for ax in axes_triplet:
        ax.set_ylim(0, max(use_max * 1.20, 0.05))
    return y_max


def plot_compare_page(genes_page, persample_by_subtype, page_idx, n_pages, out_pdf_page):
    """One page: 6 gene rows × 3 subtype columns."""
    n_genes = len(genes_page)
    n_rows = max(n_genes, 1)
    fig, axes = plt.subplots(n_rows, 3,
                              figsize=(3 * 3.0, n_rows * 1.5),
                              squeeze=False)
    for r, gene in enumerate(genes_page):
        triplet = [axes[r, 0], axes[r, 1], axes[r, 2]]
        # Use shared y-max from data in this row
        render_gene_row(triplet, gene, persample_by_subtype)
        # Gene title spanning the row — use suptitle-like position
        # We'll instead set title on the middle panel
        triplet[1].set_title(disp(gene), fontsize=11, fontweight='bold')
        # Cell type subtitles on the leftmost above bar plots? Use figure annotation per column header on top row
        # Y label on left panel
        triplet[0].set_ylabel('log-norm expr (per-mouse)', fontsize=6)

    # Add cell-type column headers above the top row only
    if n_genes > 0:
        for c, label in enumerate(SUBTYPE_LABELS):
            ax = axes[0, c]
            # Add a small text on top of axis frame
            x_center = (ax.get_xlim()[0] + ax.get_xlim()[1]) / 2
            ax.text(x_center, ax.get_ylim()[1] * 1.30, label,
                    fontsize=10, fontweight='bold', ha='center', va='bottom',
                    color={'Excitatory':'#1f77b4','Inhibitory':'#9467bd','Astrocyte':'#d62728'}[label])

    # Hide unused rows
    for r in range(n_genes, n_rows):
        for c in range(3):
            axes[r, c].axis('off')

    legend_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor='black',
                            alpha=0.8, label=g) for g in ALL_GROUPS]
    legend_handles.append(plt.Line2D([0],[0], marker='o', color='w',
                                      markerfacecolor='black', markersize=5,
                                      label='per-mouse mean'))
    legend_handles.append(plt.Line2D([0],[0], color='black', lw=0.8, label='SEM'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=5, frameon=False, fontsize=8)
    fig.suptitle(f'Excitatory | Inhibitory | Astrocyte — page {page_idx + 1} of {n_pages}',
                 fontsize=10, y=1.04)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
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
    ap.add_argument('--per-page', type=int, default=6, help='Gene rows per page')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    log_path = args.out_dir / f'{args.label}_3subtype_compare_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== 3-subtype side-by-side bar plots ===')
    logging.info('Genes listed: %d', len(GENE_LIST))

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

    is_excit = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Excitatory'] > 0))
    is_inhib = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Inhibitory'] > 0))
    is_astro = (adata.obs['celltype'].astype(str) == 'Astrocyte')

    available = [g for g in GENE_LIST if g in adata.var_names]
    missing = [g for g in GENE_LIST if g not in adata.var_names]
    logging.info('Genes plottable: %d / %d', len(available), len(GENE_LIST))
    if missing:
        logging.warning('Missing from panel: %s', missing)
    user_genes = set(USER_GENES)
    user_present = [g for g in USER_GENES if g in available]
    user_missing = [g for g in USER_GENES if g not in available]
    logging.info("User's 95 custom genes: %d present / %d missing", len(user_present), len(user_missing))
    if user_missing:
        logging.warning('User genes missing from panel: %s', user_missing)

    # Compute per-mouse means for each subtype
    persample_by_subtype = {}
    for label, mask in [('Excitatory', is_excit),
                        ('Inhibitory', is_inhib),
                        ('Astrocyte', is_astro)]:
        sub = adata[mask].copy()
        sub = sub[sub.obs['group'].astype(str).isin(ALL_GROUPS)].copy()
        sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str),
                                          categories=ALL_GROUPS, ordered=True)
        n_mice = {g: int(sub.obs.loc[sub.obs['group'] == g, 'sample_id'].nunique()) for g in ALL_GROUPS}
        logging.info('%s: %d cells; mice per group: %s', label, sub.shape[0], n_mice)
        g_arr = sub.obs['group'].astype(str).values
        s_arr = sub.obs['sample_id'].astype(str).values
        rows = []
        for gene in available:
            X = sub[:, gene].X
            v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
            df = pd.DataFrame({'expr': v, 'group': g_arr, 'sample_id': s_arr})
            mean = df.groupby(['group', 'sample_id'])['expr'].mean().reset_index()
            mean.rename(columns={'expr': 'mean_expr_log'}, inplace=True)
            mean['gene'] = gene
            rows.append(mean)
        sample_df = pd.concat(rows, ignore_index=True)
        sample_df.to_csv(
            args.out_dir / f'{args.label}_3subtype_{label}_persample_means_{args.date}.csv',
            index=False)
        persample_by_subtype[label] = sample_df

    # Plot
    per_page = args.per_page
    n_pages = (len(available) + per_page - 1) // per_page
    pdf_out = args.out_dir / f'{args.label}_3subtype_compare_BarPlots_{args.date}.pdf'
    with PdfPages(pdf_out) as pdf:
        for page_idx in range(n_pages):
            chunk = available[page_idx * per_page:(page_idx + 1) * per_page]
            plot_compare_page(chunk, persample_by_subtype, page_idx, n_pages, pdf)
    logging.info('Wrote %s (%d pages, %d genes plotted)', pdf_out, n_pages, len(available))


if __name__ == '__main__':
    main()
