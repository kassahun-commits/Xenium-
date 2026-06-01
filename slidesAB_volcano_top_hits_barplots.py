#!/usr/bin/env python3
"""
For each volcano-plot contrast, take the top 30 UP and top 30 DOWN genes
from each cell-type panel (Excitatory / Inhibitory / Astrocyte), union the
genes, and plot them as grouped-cell-type bar plots (the layout from the
prior session: each gene = one panel with 21 bars = 7 groups × 3 cell types).

Genes are organized by FUNCTIONAL CATEGORY (same 26 categories as before)
with a fall-back "Uncategorized" group for genes not in any category.

One PDF per contrast, all written into a new sub-folder `volcano_top_hits/`.

Inputs:
  * Same slide / annotation paths as before
  * Combined DE table from slidesAB_volcano_grid.py
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
    # APP added
    'APP 1':          ('APP',        1),
    'APP 2':          ('APP',        2),
    'APP 3':          ('APP',        3),
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
              'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE', 'APP']
SUBTYPES = ['Excitatory', 'Inhibitory', 'Astrocyte']
SUBTYPE_COLORS = {'Excitatory': '#1f77b4',
                  'Inhibitory': '#9467bd',
                  'Astrocyte':  '#d62728'}

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
    'Nr3c1': 'GR', 'Nr3c2': 'MR', 'Esr1': 'ERα', 'Esr2': 'ERβ',
    'Ar': 'AR', 'Pgr': 'PR', 'Thra': 'TRα', 'Thrb': 'TRβ',
    'Insr': 'INSR', 'Igf1r': 'IGF1R', 'Igf2r': 'IGF2R',
    'Lepr': 'LEPR', 'Ghr': 'GHR',
    'Crhr1': 'CRHR1', 'Crhr2': 'CRHR2', 'Oxtr': 'OXTR',
    'Avpr1b': 'V1bR', 'Mc4r': 'MC4R',
    'Atf4': 'ATF4', 'Atf6': 'ATF6', 'Ddit3': 'CHOP', 'Hspa5': 'BiP',
    'Eif2ak3': 'PERK', 'Eif2ak2': 'PKR', 'Eif2ak4': 'GCN2',
    'Eif2s1': 'eIF2α',
    'Tsc1': 'TSC1', 'Tsc2': 'TSC2', 'Mtor': 'mTOR',
    'Pmaip1': 'NOXA', 'Bbc3': 'PUMA',
    'Ntrk1': 'TrkA', 'Ntrk2': 'TrkB', 'Ntrk3': 'TrkC',
    'Ngfr': 'p75NTR', 'Ntf3': 'NT-3', 'Ntf5': 'NT-4/5',
    'Nr4a1': 'Nur77', 'Nr4a2': 'Nurr1', 'Nr4a3': 'NOR1',
}


def disp(g: str) -> str:
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


# Same 26-category map as the grouped-celltype script
CATEGORIES = {
    'Neuronal identity (sanity)': [
        'Rbfox3','Snap25','Syn1','Map2','Slc17a7','Slc17a6','Camk2a','Camk2b',
        'Satb2','Tbr1','Neurod6','Gad1','Gad2','Slc32a1','Pvalb','Sst','Vip','Reln',
    ],
    'Astrocyte identity (sanity)': [
        'Gfap','Aqp4','Slc1a3','Aldh1l1','S100b','Aldoc','Gja1','Slc1a1','Slc6a11','Glul',
    ],
    'MCT family': ['Slc16a1','Slc16a7','Slc16a3'],
    'Glucose transporters (GLUT/SGLT)': ['Slc2a1','Slc2a3','Slc2a4','Slc2a2','Slc2a5','Slc5a2'],
    'Reporter': ['GFP'],
    'Alcohol-metabolizing enzymes': ['Adh1','Adh5','Aldh2','Aldh1a2','Cyp2e1','Cat'],
    'Acetate / Ac-CoA / lactate': ['Acss1','Acss2','Acly','Ldha','Ldhb'],
    'Glycolysis / TCA / energy': [
        'Gapdh','Pgk1','Pkm','Hk1','Hk2','Pfkm','Aldoa','Eno1','Eno2',
        'Idh1','Idh2','Cs','Pdha1','Pdhb',
    ],
    'One-carbon / SAM / methylation metabolism': [
        'Ahcy','Chdh','Gnmt','Mat1a','Mat2a','Mat2b','Mtap','Mthfd1','Mthfd2',
        'Mtr','Mtrr','Pemt','Phgdh','Psph','Shmt2','Tyms',
    ],
    'Lipid / fatty-acid metabolism': ['Acat1','Acat2','Crat','Crot','Elovl5'],
    'Solute transporters (misc)': ['Slc4a4','Slc5a8','Slc10a2','Slc13a5',
                                    'Slc25a32','Slc25a40','Bsg','Emb','Abcc3'],
    'Oxidative stress / detox': [
        'Sod1','Sod2','Gpx1','Gpx4','Nfe2l2','Nfe2l1','Keap1','Hmox1',
        'Sirt3','Sirt5','Gstp1',
    ],
    'Metabolic / ER stress / UPR / ISR': [
        'Atf4','Atf6','Ddit3','Hspa5','Eif2ak3','Eif2ak2','Eif2ak4','Eif2s1',
        'Sesn1','Sesn2','Tsc1','Tsc2','Mtor','Rheb','Foxo1','Foxo3','Foxo4',
        'Pmaip1','Bbc3','Trib3',
    ],
    'Hormone receptors': [
        'Nr3c1','Nr3c2','Esr1','Esr2','Ar','Pgr','Thra','Thrb','Insr','Igf1r',
        'Igf2r','Lepr','Ghr','Crhr1','Crhr2','Oxtr','Avpr1b','Mc4r',
    ],
    'IEGs — Fos family': ['Fos','Fosb','Fosl1','Fosl2'],
    'IEGs — Jun family': ['Jun','Junb'],
    'IEGs — Egr family': ['Egr1','Egr2','Egr3','Egr4'],
    'IEGs — Nr4a family': ['Nr4a1','Nr4a2','Nr4a3'],
    'IEGs — activity-regulated': [
        'Arc','Bdnf','Homer1','Npas4','Per1','Per2','Dusp1','Dusp6',
        'Plk2','Plk3','Nptx1','Nptx2','Trib2',
    ],
    'IEGs — other / activity-related': [
        'Atf3','Btg2','Gadd45b','Bhlhe40','Klf10','Ccn1','Rasgrf1',
    ],
    'Neurotrophic factors': ['Ntf3','Ntf5','Ngf','Gdnf','Cntf','Lif','Mst1','Mstn','Igf1','Igf2'],
    'Neurotrophic receptors': ['Ntrk1','Ntrk2','Ntrk3','Ngfr','Gfra1','Gfra2','Gfra3','Ret','Lifr'],
    'Chromatin / histone modifiers': [
        'Brd9','Chd3','Chd8','Dpf1','Ino80','Kat7','Kmt5a','Mbd1','Mbd4',
        'Nsd2','Nsd3','Prdm2','Prdm8','Prmt2','Prmt6','Prmt8','Prmt9',
        'Setdb2','Ss18','Suv39h2',
    ],
    'Synaptic / structural': [
        'Grik5','Rgs14','Pcp4','Cacng5','Cacnb4','Dlg2','Lrrtm2','Lrrtm4',
        'Pclo','Unc13a','Tagln3',
    ],
    'RNA-binding / nuclear / splicing': [
        'Celf2','Ccnt2','Cdk2ap1','Hmgb3','Hmgn3','Ilf3','Kcnh7','Malat1',
        'Mbnl2','Pnisr','Prpf39','Rbm34','Rsrc2','Tia1','Zranb2',
    ],
    'Other / signaling': ['Cfhr1','Abhd14b','Dgkb','Gng12'],
}

COMPARISONS = [
    ('EtOH_veh',    'H2O_veh',  'Acute_alcohol__EtOHveh_vs_H2Oveh',
     'Acute alcohol — EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic_alcohol__ChronicEtOH_vs_H2Oveh',
     'Chronic alcohol — ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'DrugOnly__H2OMCT1i_vs_H2Oveh',
     'Drug-only — H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'AlcoholDrug__EtOHMCT1i_vs_H2Oveh',
     'Alcohol + drug — EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'MCT1iRescue__EtOHMCT1i_vs_EtOHveh',
     'MCT1i rescue — EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2AOE_vs_MAT2ACM',
     'MAT2A overexpression — MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP_vs_H2Oveh',
     'APP — APP vs H2O_veh'),
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


def plot_one_gene_panel(ax, gene, persample_by_subtype, direction_tag=None):
    """One panel with 21 bars: 7 group clusters × 3 cell-types per cluster."""
    group_positions = np.arange(len(ALL_GROUPS))
    offsets = {'Excitatory': -0.27, 'Inhibitory': 0.0, 'Astrocyte': 0.27}
    bar_width = 0.25
    y_max = 0.0
    for gi, grp in enumerate(ALL_GROUPS):
        base_x = group_positions[gi]
        for subtype in SUBTYPES:
            d = persample_by_subtype[subtype]
            d = d[(d['gene'] == gene) & (d['group'] == grp)]
            if len(d) == 0: continue
            means = d['mean_expr_log'].values
            mh = means.mean()
            sem = means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1 else 0.0
            x = base_x + offsets[subtype]
            ax.bar(x, mh, width=bar_width, color=SUBTYPE_COLORS[subtype],
                   edgecolor='black', linewidth=0.4, alpha=0.85, zorder=2)
            if sem > 0:
                ax.errorbar(x, mh, yerr=sem, color='black', capsize=1.5, lw=0.6, zorder=3)
            np.random.seed(hash((gene, subtype, grp)) % (2**32))
            jit = np.random.uniform(-0.06, 0.06, size=len(means))
            ax.scatter([x] * len(means) + jit, means,
                       c='black', s=5, zorder=4, edgecolor='white', linewidth=0.3)
            y_max = max(y_max, mh + sem, float(means.max()) if len(means) else 0)
    for sep in group_positions[:-1] + 0.5:
        ax.axvline(sep, color='lightgray', lw=0.5, alpha=0.6, zorder=1)
    ax.set_xticks(group_positions)
    ax.set_xticklabels(ALL_GROUPS, rotation=45, ha='right', fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_xlim(-0.55, len(ALL_GROUPS) - 0.45)
    ax.set_ylim(0, max(y_max * 1.20, 0.05))
    title = disp(gene)
    if direction_tag:
        title += f'   [{direction_tag}]'
    ax.set_title(title, fontsize=10, fontweight='bold')


def make_contrast_pdf(test, ref, slug, pretty_label,
                       de_long, persample, all_panel_genes, out_pdf,
                       n_top=30, padj_thresh=1e-3):
    """One PDF per contrast: top 30 up + top 30 down genes per cell type, unioned, functionally grouped."""
    # Subset DE to this contrast
    contrast_de = de_long[(de_long['test'] == test) & (de_long['reference'] == ref)]
    if contrast_de.empty:
        logging.warning('No DE rows for %s vs %s', test, ref); return None

    # Top hits per cell type
    per_ct_top = {}
    for ct in SUBTYPES:
        d = contrast_de[(contrast_de['subtype'] == ct) & (contrast_de['padj'] < padj_thresh)]
        up = d.nlargest(n_top, 'logfc')
        dn = d.nsmallest(n_top, 'logfc')
        per_ct_top[ct] = {'up': up, 'down': dn}

    # Direction map: each gene tagged with which cell types had it up/down
    direction = {}  # gene -> set of ('E/I/A: up' or 'E/I/A: dn')
    ct_short = {'Excitatory': 'E', 'Inhibitory': 'I', 'Astrocyte': 'A'}
    for ct in SUBTYPES:
        for g in per_ct_top[ct]['up']['gene'].tolist():
            direction.setdefault(g, set()).add(f"{ct_short[ct]}:↑")
        for g in per_ct_top[ct]['down']['gene'].tolist():
            direction.setdefault(g, set()).add(f"{ct_short[ct]}:↓")
    all_top = sorted(direction.keys())
    n_total = len(all_top)
    logging.info('  %s: %d unique top hits (top %d up/down per cell type)', slug, n_total, n_top)

    # Functional categorization
    seen = set()
    cat_genes = {}
    for cat, gs in CATEGORIES.items():
        keep = [g for g in gs if g in direction and g not in seen]
        seen.update(keep)
        if keep:
            cat_genes[cat] = keep
    uncategorized = [g for g in all_top if g not in seen]
    if uncategorized:
        cat_genes['Uncategorized'] = uncategorized

    # Also save the top hits table
    rows = []
    for ct in SUBTYPES:
        for direction_label, dframe in [('up', per_ct_top[ct]['up']),
                                          ('down', per_ct_top[ct]['down'])]:
            for _, r in dframe.iterrows():
                rows.append({
                    'subtype': ct, 'direction': direction_label,
                    'gene': r['gene'], 'logfc': r['logfc'], 'padj': r['padj'],
                    'rank_in_subtype': None,  # filled below
                })
    top_table = pd.DataFrame(rows)

    PANELS_PER_PAGE = 6  # 2 cols × 3 rows
    page_idx = 0
    total_pages = 0
    for gs in cat_genes.values():
        total_pages += max(1, (len(gs) + PANELS_PER_PAGE - 1) // PANELS_PER_PAGE)

    with PdfPages(out_pdf) as pdf:
        # Cover page
        cover_fig = plt.figure(figsize=(11, 8.5))
        cover_fig.text(0.5, 0.93, 'Top hits — grouped-cell-type bar plots',
                       fontsize=16, ha='center', fontweight='bold')
        cover_fig.text(0.5, 0.88, pretty_label, fontsize=12, ha='center', style='italic')
        meta = (f"Top {n_top} UP + Top {n_top} DOWN genes per cell type (padj < {padj_thresh:g}, "
                f"ranked by log2 fold change), union'd across the 3 cell types.\n\n"
                f"Total unique top hits: {n_total} genes.\n\n"
                "Layout per panel: 7 condition-group clusters along x; within each cluster,\n"
                "three bars side-by-side colored by cell type — blue=Excitatory,\n"
                "purple=Inhibitory, red=Astrocyte. Bar height = mean of per-mouse means,\n"
                "error bars = SEM, dots = individual mouse means.\n\n"
                "Title tag in brackets shows which cell type(s) had this gene in the\n"
                "top hits and the direction (↑ up in test, ↓ up in reference):\n"
                "    [E:↑ A:↓]  =  up in Excitatory's volcano, down in Astrocyte's\n\n"
                "Genes are organized into 26 functional categories.\n\n"
                "Per-cell-type top-hit gene table is saved alongside as a CSV."
                "\n\nStatistics: cell-level Wilcoxon. Treat directions / log2FC as\n"
                "trustworthy; treat adj p-values as ranking statistics (pseudoreplication).")
        cover_fig.text(0.08, 0.78, meta, fontsize=10, va='top', ha='left', wrap=True)
        cover_fig.text(0.08, 0.20, f"Per-cell-type breakdown (top {n_top} each direction):",
                       fontsize=11, fontweight='bold')
        y_pos = 0.16
        for ct in SUBTYPES:
            up_n = len(per_ct_top[ct]['up'])
            dn_n = len(per_ct_top[ct]['down'])
            cover_fig.text(0.10, y_pos,
                           f"{ct:12s}  up: {up_n}   down: {dn_n}",
                           fontsize=10, fontfamily='monospace')
            y_pos -= 0.025
        pdf.savefig(cover_fig, bbox_inches='tight'); plt.close(cover_fig)

        # Plot pages, organized by category
        for cat, genes in cat_genes.items():
            n_pages = max(1, (len(genes) + PANELS_PER_PAGE - 1) // PANELS_PER_PAGE)
            for p in range(n_pages):
                chunk = genes[p * PANELS_PER_PAGE:(p + 1) * PANELS_PER_PAGE]
                fig, axes = plt.subplots(3, 2, figsize=(12, 11.5), squeeze=False)
                for i, g in enumerate(chunk):
                    r_, c_ = divmod(i, 2)
                    dtag = ' '.join(sorted(direction.get(g, [])))
                    plot_one_gene_panel(axes[r_, c_], g, persample, direction_tag=dtag)
                    axes[r_, c_].set_ylabel('log-norm expr (per-mouse)', fontsize=7)
                for j in range(len(chunk), 6):
                    r_, c_ = divmod(j, 2)
                    axes[r_, c_].axis('off')

                legend = [Patch(facecolor=SUBTYPE_COLORS[s], edgecolor='black',
                                alpha=0.85, label=s) for s in SUBTYPES]
                legend += [plt.Line2D([0], [0], marker='o', color='w',
                                       markerfacecolor='black', markersize=6,
                                       label='per-mouse mean'),
                           plt.Line2D([0], [0], color='black', lw=0.8, label='SEM')]
                fig.legend(handles=legend, loc='upper center',
                           bbox_to_anchor=(0.5, 1.005), ncol=5, frameon=False, fontsize=9)
                page_idx += 1
                fig.suptitle(f"{pretty_label}\n"
                             f"Category: {cat}  (page {page_idx} of {total_pages}; "
                             f"cat sub-page {p + 1}/{n_pages})",
                             fontsize=11, y=1.02)
                plt.tight_layout(rect=[0, 0, 1, 0.98])
                pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

    return top_table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-a-ann', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--de-csv', required=True, type=Path,
                    help='Combined DE CSV from slidesAB_volcano_grid.py')
    ap.add_argument('--out-dir', required=True, type=Path,
                    help='Parent output directory (creates volcano_top_hits/ inside)')
    ap.add_argument('--label', default='SlidesAB')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--n-top', type=int, default=30)
    ap.add_argument('--padj-thresh', type=float, default=1e-3)
    args = ap.parse_args()

    # Make new sub-folder
    vp_folder = args.out_dir / 'volcano_top_hits'
    vp_folder.mkdir(parents=True, exist_ok=True)

    log_path = vp_folder / f'{args.label}_volcano_top_hits_log.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s',
                        handlers=[logging.FileHandler(log_path, mode='w'),
                                  logging.StreamHandler(sys.stdout)])
    logging.info('=== Volcano top-hits bar plots ===')

    # Load combined slides + per-sample means (per cell type)
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

    # Load DE
    de_long = pd.read_csv(args.de_csv)
    logging.info('Loaded DE: %d rows', len(de_long))

    # Collect all genes that may appear in any top hit
    needed_genes = set()
    for test, ref, _, _ in COMPARISONS:
        cde = de_long[(de_long['test'] == test) & (de_long['reference'] == ref)]
        for ct in SUBTYPES:
            d = cde[(cde['subtype'] == ct) & (cde['padj'] < args.padj_thresh)]
            needed_genes.update(d.nlargest(args.n_top, 'logfc')['gene'].tolist())
            needed_genes.update(d.nsmallest(args.n_top, 'logfc')['gene'].tolist())
    available = [g for g in needed_genes if g in adata.var_names]
    logging.info('Unique genes across all contrasts top hits: %d (in panel: %d)',
                 len(needed_genes), len(available))

    # Compute per-mouse means once per subtype, only for needed genes
    persample = {}
    for subtype, mask in [('Excitatory', is_excit),
                          ('Inhibitory', is_inhib),
                          ('Astrocyte',  is_astro)]:
        sub = adata[mask].copy()
        sub = sub[sub.obs['group'].astype(str).isin(ALL_GROUPS)].copy()
        sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str),
                                          categories=ALL_GROUPS, ordered=True)
        n_mice = {g: int(sub.obs.loc[sub.obs['group'] == g, 'sample_id'].nunique()) for g in ALL_GROUPS}
        logging.info('%s: %d cells; mice/group: %s', subtype, sub.shape[0], n_mice)
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
        persample[subtype] = pd.concat(rows, ignore_index=True)

    # Per-contrast PDFs
    top_tables = []
    for test, ref, slug, pretty in COMPARISONS:
        out_pdf = vp_folder / f'{args.label}_TopHits_{slug}_{args.date}.pdf'
        top_table = make_contrast_pdf(test, ref, slug, pretty,
                                       de_long, persample, available, out_pdf,
                                       n_top=args.n_top, padj_thresh=args.padj_thresh)
        if top_table is not None:
            top_table['contrast'] = slug
            top_table['test'] = test; top_table['reference'] = ref
            top_tables.append(top_table)
        logging.info('Wrote %s', out_pdf)

    if top_tables:
        all_top = pd.concat(top_tables, ignore_index=True)
        all_top.to_csv(vp_folder / f'{args.label}_TopHits_AllContrasts_{args.date}.csv',
                       index=False)
        logging.info('Wrote per-cell-type top-hit table CSV')


if __name__ == '__main__':
    main()
