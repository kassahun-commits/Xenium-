#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled — FOLD-CHANGE curated dot plot, 3 cell types.

New version of the curated dot plot (cf. FINAL DOT BLOT
SlidesAB_CuratedDotPlot_NeuronAstrocyte) with these changes:

  * CONSTANT dot size (no encoding of fraction expressing).
  * COLOR encodes DIRECTION of change: log2 fold change of each treatment
    group vs the H2O_veh vehicle control (diverging blue-white-red, centred 0).
  * THREE cell-type panels side by side: Excitatory / Inhibitory / Astrocyte.
  * Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).
  * Genes grouped into the 17 curated functional categories, PLUS a final
    section with the user-provided additional gene list (97 genes + GFP = "98").

Columns of every panel (all vs H2O_veh):
    MCT1i        = H2O_MCT1i   vs H2O_veh
    EtOH         = EtOH_veh    vs H2O_veh
    EtOH+MCT1i   = EtOH_MCT1i  vs H2O_veh
    Chronic EtOH = ChronicEtOH vs H2O_veh

Fold change is computed per cell type as a pseudobulk ratio:
    log2FC = log2( (mean_norm_expr[group] + 1) / (mean_norm_expr[H2O_veh] + 1) )
where mean_norm_expr is the mean of linear (expm1 of log-normalised) expression
across all cells of that cell type in that group (Slides A + B pooled).
A group x cell-type cell-count guard (>= MIN_CELLS) blanks under-powered cells
(notably Inhibitory in V2) instead of plotting noise.

Cell typing follows the master build: argmax marker-module score for
Neuron/Astrocyte/Oligodendrocyte/Microglia, then Excit/Inhib split within
neurons (excitatory precedence if both score > 0).

All paths via CLI (MEWS Lab rule). Vector PDFs, editable text (pdf.fonttype=42).
One source CSV per version is written alongside the figures.
"""
from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import anndata as ad
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ----------------------------- parameters ----------------------------------
# Columns (left -> right): (id, test group, reference group, two-line label).
# Columns 1-4 are vs the H2O_veh vehicle control; the final column is the direct
# MCT1i-rescue contrast EtOH_MCT1i vs EtOH_veh (does blocking monocarboxylate
# uptake reverse the ethanol effect?).
COLS = [
    ('MCT1i_vs_H2O',       'H2O_MCT1i',   'H2O_veh',  'MCT1i\nvs H2O'),
    ('EtOH_vs_H2O',        'EtOH_veh',    'H2O_veh',  'EtOH\nvs H2O'),
    ('EtOHMCT1i_vs_H2O',   'EtOH_MCT1i',  'H2O_veh',  'EtOH+MCT1i\nvs H2O'),
    ('ChronicEtOH_vs_H2O', 'ChronicEtOH', 'H2O_veh',  'Chronic EtOH\nvs H2O'),
    ('EtOHMCT1i_vs_EtOH',  'EtOH_MCT1i',  'EtOH_veh', 'EtOH+MCT1i\nvs EtOH'),
]
ALCOHOL_GROUPS = sorted({g for _id, t, r, _l in COLS for g in (t, r)})
PANELS = ['Excitatory', 'Inhibitory', 'Astrocyte']     # 'subtype' column
CLIP = 2.0          # log2FC colour clip (+/-)
DOT_SIZE = 60.0     # constant marker size
MIN_CELLS = 10      # min cells per (group x cell type) to plot a value
PSEUDO = 1.0        # pseudocount on linear normalised scale

# ----------------------------- gene sets ------------------------------------
GENE_CATEGORIES = {
    'Cell-type identity':
        ['Rbfox3', 'Snap25', 'Syn1', 'Map2', 'Gfap', 'Aqp4', 'Slc1a3', 'S100b', 'Aldh1l1'],
    'Alcohol & acetaldehyde metab':
        ['Adh1', 'Adh5', 'Aldh2', 'Aldh1a1', 'Cyp2e1', 'Cat', 'Akr1a1'],
    'Acetate / Ac-CoA metab':
        ['Acss1', 'Acss2', 'Acly', 'Acaca', 'Acacb'],
    'Lactate / MCT axis':
        ['Slc16a1', 'Slc16a3', 'Slc16a7', 'Ldha', 'Ldhb'],
    'Glucose transport':
        ['Slc2a1', 'Slc2a2', 'Slc2a3', 'Slc2a4', 'Slc2a5', 'Slc5a2'],
    'Glycolysis / TCA':
        ['Gapdh', 'Pgk1', 'Pkm', 'Hk1', 'Hk2', 'Pfkm', 'Aldoa', 'Eno1', 'Eno2', 'Idh1', 'Idh2', 'Cs'],
    'SAM / MAT2A axis':
        ['Mat1a', 'Mat2a', 'Mat2b', 'Ahcy', 'Mthfr', 'Cbs'],
    'Glutathione / oxidative stress':
        ['Gclc', 'Gclm', 'Gpx1', 'Gpx4', 'Gss', 'Gstm1', 'Sod1', 'Sod2', 'Hmox1', 'Nfe2l2'],
    'Stress / UPR':
        ['Atf3', 'Atf4', 'Ddit3', 'Hspa1a', 'Hsp90aa1', 'Hspb1', 'Trp53'],
    'Activity / IEGs':
        ['Fos', 'Jun', 'Junb', 'Arc', 'Egr1', 'Egr2', 'Npas4', 'Homer1', 'Bdnf'],
    'Glutamate signaling':
        ['Slc17a7', 'Slc17a6', 'Slc1a2', 'Grin1', 'Grin2a', 'Grin2b', 'Gria1', 'Gria2',
         'Grm1', 'Grm2', 'Grm5', 'Camk2a', 'Camk2b'],
    'GABA signaling':
        ['Gad1', 'Gad2', 'Gabra1', 'Gabra2', 'Gabbr1', 'Slc32a1'],
    'Endocannabinoid':
        ['Cnr1', 'Faah', 'Mgll', 'Daglb'],
    'Neuroinflammation':
        ['Tnf', 'Il1b', 'Il6', 'Il10', 'Tlr4', 'Nfkb1', 'Stat3'],
    'Glutamate transporters (EAAT/Slc1a)':
        ['Slc1a1', 'Slc1a3'],
    'HATs (lysine acetyltransferases)':
        ['Ep300', 'Crebbp', 'Kat2a', 'Kat2b', 'Kat5', 'Kat6a', 'Kat6b', 'Kat7', 'Kat8',
         'Clock', 'Ncoa1', 'Ncoa2', 'Ncoa3'],
    'HDACs / Sirtuins':
        ['Hdac1', 'Hdac2', 'Hdac3', 'Hdac4', 'Hdac5', 'Hdac6', 'Hdac8', 'Hdac9',
         'Sirt1', 'Sirt2', 'Sirt3', 'Sirt4', 'Sirt5', 'Sirt6', 'Sirt7'],
}

# Additional gene list provided by the user (97 pasted + GFP = "98").
EXTRA_SECTION = 'Additional user gene list (n=98)'
EXTRA_GENES = [
    'GFP',
    'Ahcy', 'Atf3', 'Brd9', 'Btg2', 'Chd3', 'Chd8', 'Chdh', 'Dpf1', 'Fos', 'Fosb',
    'Gadd45b', 'Gnmt', 'Ino80', 'Jun', 'Junb', 'Kat7', 'Kmt5a', 'Mat1a', 'Mat2a',
    'Mat2b', 'Mbd1', 'Mbd4', 'Mtap', 'Mthfd1', 'Mthfd2', 'Mtr', 'Mtrr', 'Nsd2', 'Nsd3',
    'Pemt', 'Phgdh', 'Prdm2', 'Prdm8', 'Prmt2', 'Prmt6', 'Prmt8', 'Prmt9', 'Psph',
    'Setdb2', 'Shmt2', 'Slc25a32', 'Ss18', 'Suv39h2', 'Tyms', 'Grik5', 'Rgs14', 'Pcp4',
    'Cacng5', 'Cacnb4', 'Ccnt2', 'Cdk2ap1', 'Celf2', 'Dlg2', 'Hmgb3', 'Hmgn3', 'Ilf3',
    'Kcnh7', 'Lrrtm2', 'Lrrtm4', 'Malat1', 'Mbnl2', 'Pclo', 'Pnisr', 'Prpf39', 'Rbm34',
    'Rsrc2', 'Slc25a40', 'Tia1', 'Zranb2', 'Slc5a8', 'Emb', 'Bsg', 'Slc4a4', 'Crat',
    'Crot', 'Slc10a2', 'Acat1', 'Acat2', 'Slc16a7', 'Gapdh', 'Slc13a5', 'Pdha1', 'Pdhb',
    'Gpx4', 'Gpx1', 'Glul', 'Adh1', 'Unc13a', 'Abcc3', 'Abhd14b', 'Tagln3', 'Nptx1',
    'Pfkm', 'Cfhr1', 'Dgkb', 'Gng12', 'Elovl5',
]

ALIAS = {'Slc2a1': 'GLUT1', 'Slc2a3': 'GLUT3', 'Cat': 'catalase'}

# ----------------------------- loaders (master build) ------------------------
LABEL_MAP = {
    'Chronic EtOH 1': ('ChronicEtOH', 1), 'MCT1i+EtOH 1': ('EtOH_MCT1i', 1),
    'MCT1i+Veh 1': ('H2O_MCT1i', 1), 'MAT2A CM 1': ('MAT2A_CM', 1),
    'MAT2A CM 2': ('MAT2A_CM', 2), 'MAT2A CM 3': ('MAT2A_CM', 3),
    'MAT2A OE 1': ('MAT2A_OE', 1), 'MAT2A OE 2': ('MAT2A_OE', 2),
    'MAT2A OE 3': ('MAT2A_OE', 3), 'Veh+EtOH 1': ('EtOH_veh', 1),
    'veh +EtOH 2': ('EtOH_veh', 2), 'veh+H20 1': ('H2O_veh', 1),
    'EtOH +MCT1i 2': ('EtOH_MCT1i', 2), 'EtOH +MCT1i 3': ('EtOH_MCT1i', 3),
    'H20+veh 2': ('H2O_veh', 2), 'H20 +MCT1i 2': ('H2O_MCT1i', 2),
    'Chronic EtOH 2': ('ChronicEtOH', 2), 'chronic EtOH 3': ('ChronicEtOH', 3),
    'EtOH+MCT1i 1': ('EtOH_MCT1i', 1), 'H20+MCT1i 1': ('H2O_MCT1i', 1),
    'H20 +MCT1i 1': ('H2O_MCT1i', 1), 'H20 +veh 1': ('H2O_veh', 1),
    'EtOH +veh 1': ('EtOH_veh', 1), 'APP 1': ('APP', 1),
    'APP 2': ('APP', 2), 'APP 3': ('APP', 3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']
INHIB_MARKERS = ['Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln', 'Lhx6']


def load_polygons(csv):
    df = pd.read_csv(csv, comment='#')
    return {name: sub[['X', 'Y']].values.astype(float)
            for name, sub in df.groupby('Selection', sort=False)}


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


def cell_type_data(adata):
    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gs_present, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    scores = adata.obs[[f'score_{ct}' for ct in MARKERS]]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_ct[scores.max(axis=1) <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best_ct, categories=list(MARKERS.keys()) + ['Unclassified'])
    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory', random_state=0, n_bins=25)
    inhib_present = [g for g in INHIB_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=inhib_present, score_name='score_Inhibitory', random_state=0, n_bins=25)
    is_neuron = adata.obs['celltype'].astype(str) == 'Neuron'
    is_excit = is_neuron & (adata.obs['score_Excitatory'] > 0)
    is_inhib = is_neuron & (adata.obs['score_Inhibitory'] > 0)
    is_astro = adata.obs['celltype'].astype(str) == 'Astrocyte'
    subtype = pd.Series('Other', index=adata.obs.index)
    subtype[is_inhib] = 'Inhibitory'
    subtype[is_excit] = 'Excitatory'  # excit precedence
    subtype[is_astro] = 'Astrocyte'
    adata.obs['subtype'] = pd.Categorical(subtype, categories=['Excitatory', 'Inhibitory', 'Astrocyte', 'Other'])
    return adata


def process_version(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    ada_a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    ada_b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return cell_type_data(adata)


def dense_col(adata, gene):
    X = adata[:, gene].X
    return np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()


def disp(g):
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


# ----------------------------- row layout -----------------------------------
def build_rows():
    """Ordered (section, gene) rows. Functional categories deduped across
    categories (keep first occurrence). Extra section kept whole (deduped only
    within itself), so it can intentionally repeat earlier genes."""
    rows = []
    seen = set()
    for section, genes in GENE_CATEGORIES.items():
        for g in genes:
            if g in seen:
                continue
            seen.add(g)
            rows.append((section, g))
    seen_extra = set()
    for g in EXTRA_GENES:
        if g in seen_extra:
            continue
        seen_extra.add(g)
        rows.append((EXTRA_SECTION, g))
    return rows


# ----------------------------- fold-change ----------------------------------
def fc_table(adata, rows):
    """Return long-form DataFrame of log2FC per (gene, cell type, column) where
    each column carries its own (test, reference) group pair; plus n_cells /
    mean expr, for the source CSV and plotting."""
    sub = adata[(adata.obs['subtype'].astype(str).isin(PANELS)) &
                (adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS))].copy()
    var = set(sub.var_names)
    genes_needed = sorted({g for _, g in rows})
    present = [g for g in genes_needed if g in var]

    # mean linear normalised expr and cell counts per (gene, ct, group)
    grp = sub.obs['group'].astype(str).values
    ctv = sub.obs['subtype'].astype(str).values
    n_cells = (pd.DataFrame({'ct': ctv, 'group': grp})
               .groupby(['ct', 'group']).size())

    mean_lin = {}    # (gene, ct, group) -> mean linear normalised expr
    for gene in present:
        lin = np.expm1(dense_col(sub, gene))
        d = pd.DataFrame({'lin': lin, 'ct': ctv, 'group': grp})
        m = d.groupby(['ct', 'group'])['lin'].mean()
        for (ct, group), val in m.items():
            mean_lin[(gene, ct, group)] = val

    out = []
    for section, gene in rows:
        in_panel = gene in var
        for ct in PANELS:
            for cid, test, ref, _lab in COLS:
                n_t = int(n_cells.get((ct, test), 0))
                n_r = int(n_cells.get((ct, ref), 0))
                m_t = mean_lin.get((gene, ct, test), np.nan) if in_panel else np.nan
                m_r = mean_lin.get((gene, ct, ref), np.nan) if in_panel else np.nan
                if (not in_panel or n_t < MIN_CELLS or n_r < MIN_CELLS or
                        pd.isna(m_t) or pd.isna(m_r)):
                    l2fc = np.nan
                else:
                    l2fc = float(np.log2((m_t + PSEUDO) / (m_r + PSEUDO)))
                out.append({'section': section, 'gene': gene, 'in_panel': in_panel,
                            'cell_type': ct, 'column': cid,
                            'test_group': test, 'ref_group': ref,
                            'n_test': n_t, 'n_ref': n_r,
                            'mean_norm_expr_test': m_t, 'mean_norm_expr_ref': m_r,
                            'log2fc': l2fc})
    return pd.DataFrame(out)


# ----------------------------- drawing --------------------------------------
def draw_page(pdf, adata, version, rows, fc_df):
    n_rows = len(rows)
    y_of = {(sec, g): i for i, (sec, g) in enumerate(rows)}  # 0 = top

    # section boundaries (row index where each section starts) + centres
    sec_order, sec_start = [], {}
    for i, (sec, _g) in enumerate(rows):
        if sec not in sec_start:
            sec_start[sec] = i
            sec_order.append(sec)
    sec_bounds = [sec_start[s] for s in sec_order] + [n_rows]

    ncol = len(COLS)
    panel_counts = (adata.obs[adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS)]
                    ['subtype'].astype(str).value_counts().to_dict())

    fig_w = 3.6 + len(PANELS) * (0.62 * ncol + 0.8)
    fig_h = max(8.0, 0.135 * n_rows + 2.4)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(fig_w, fig_h), sharey=True,
                             gridspec_kw={'wspace': 0.12})
    if len(PANELS) == 1:
        axes = [axes]

    norm = plt.Normalize(vmin=-CLIP, vmax=CLIP)
    cmap = plt.cm.RdBu_r

    for ax, ct in zip(axes, PANELS):
        d = fc_df[fc_df['cell_type'] == ct]
        # quick lookup keyed on the column id
        lut = {(r.gene, r.section, r.column): r.log2fc for r in d.itertuples(index=False)}
        for (sec, gene), yi in y_of.items():
            for j, (cid, _t, _r, _lab) in enumerate(COLS):
                val = lut.get((gene, sec, cid), np.nan)
                if pd.isna(val):
                    continue
                ax.scatter(j, yi, s=DOT_SIZE, c=[val], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.3, zorder=3)
        ax.set_xlim(-0.6, ncol - 0.4)
        ax.set_ylim(n_rows - 0.5, -0.5)   # row 0 at top
        col_labels = [lab for _c, _t, _r, lab in COLS]
        ax.set_xticks(range(ncol))
        ax.set_xticklabels(col_labels, rotation=0, ha='center', va='top', fontsize=6.0)
        # repeat column labels at the top of the (very tall) panel for readability
        secax = ax.secondary_xaxis('top')
        secax.set_xticks(range(ncol))
        secax.set_xticklabels(col_labels, rotation=0, ha='center', va='bottom', fontsize=6.0)
        secax.tick_params(length=0)
        for sp in secax.spines.values():
            sp.set_visible(False)
        # mark the rescue column (different reference) with a dashed divider
        ax.axvline(ncol - 1.5, color='0.4', lw=0.7, ls='--', zorder=1)
        ax.set_title(f'{ct}\n(n={panel_counts.get(ct, 0):,} cells)',
                     fontsize=9, fontweight='bold', pad=24)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # section separators
        for b in sec_bounds[1:-1]:
            ax.axhline(b - 0.5, color='0.6', lw=0.6, zorder=1)
        ax.set_axisbelow(True)
        ax.grid(axis='x', color='0.92', lw=0.4, zorder=0)

    # y labels (gene names) on leftmost axis; mark genes absent from panel
    present_genes = set(fc_df.loc[fc_df['in_panel'], 'gene'])
    ylabels, ycolors = [], []
    for _s, g in rows:
        if g in present_genes:
            ylabels.append(disp(g)); ycolors.append('black')
        else:
            ylabels.append(disp(g) + ' †'); ycolors.append('0.6')
    fs = max(3.0, min(6.0, 560 / max(n_rows, 1)))
    axes[0].set_yticks(range(n_rows))
    axes[0].set_yticklabels(ylabels, fontsize=fs)
    for tick, col in zip(axes[0].get_yticklabels(), ycolors):
        tick.set_color(col)

    # section labels on the right of the rightmost panel
    rax = axes[-1]
    for s, e, sec in zip(sec_bounds[:-1], sec_bounds[1:], sec_order):
        yc = (s + e - 1) / 2.0
        rax.text(ncol - 0.30, yc, sec, fontsize=6.2, va='center', ha='left',
                 rotation=0, clip_on=False, color='0.15')

    # colour bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.16)
    cbar.set_label('log2 fold change (per-column reference)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    n_absent = sum(1 for _s, g in rows if g not in present_genes)
    fig.suptitle(
        f'Slides A+B pooled — fold-change curated dot plot — {version}\n'
        f'colour = log2FC (clipped ±{CLIP:g}); cols 1-4 vs H2O_veh, last col = '
        f'EtOH+MCT1i vs EtOH (rescue, dashed); constant dot size; 3 cell types\n'
        f'(blank cell = gene absent/under-powered; † = not in panel, n={n_absent}; '
        f'min {MIN_CELLS} cells per group×cell-type)',
        fontsize=8.5, y=1.0)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def cover(title, body, path):
    d = fitz.open(); p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 54, 738, 150), title, fontsize=16,
                     fontname='Helvetica-Bold', align=1)
    p.insert_textbox(fitz.Rect(54, 168, 738, 560), body, fontsize=11,
                     fontname='Helvetica')
    d.save(str(path)); d.close()


# ----------------------------- main -----------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-b-ann', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--clip', type=float, default=2.0,
                    help='log2FC colour clip (+/-); default 2.0. Non-default '
                         'values add a _clip<val> tag to output filenames.')
    args = ap.parse_args()

    global CLIP
    CLIP = args.clip
    clip_tag = '' if abs(args.clip - 2.0) < 1e-9 else f"_clip{('%g' % args.clip).replace('.', 'p')}"

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)
    logging.info('log2FC colour clip = +/-%g (filename tag: %r)', CLIP, clip_tag)

    rows = build_rows()
    logging.info('Total rows: %d (%d functional + %d extra-section)',
                 len(rows), sum(1 for s, _ in rows if s != EXTRA_SECTION),
                 sum(1 for s, _ in rows if s == EXTRA_SECTION))

    tmp = args.outdir / '_tmp_fc_dot'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    versions = {
        'V1': (args.v1_slide_a_ann, args.v1_slide_b_ann),
        'V2': (args.v2_slide_a_ann, args.v2_slide_b_ann),
    }
    with PdfPages(pages_pdf) as pdf:
        for version, (a_ann, b_ann) in versions.items():
            logging.info('=== Processing %s (pool A+B) ===', version)
            adata = process_version(args.slide_a_dir, a_ann, args.slide_b_dir, b_ann)
            comp = (adata.obs[adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS)]
                    ['subtype'].astype(str).value_counts().to_dict())
            logging.info('%s subtype composition (alcohol groups): %s', version, comp)
            fc_df = fc_table(adata, rows)
            csv_path = args.outdir / f'SlidesAB_FCcuratedDotplot_3CT_{version}{clip_tag}_{today}.csv'
            fc_df.to_csv(csv_path, index=False)
            logging.info('Wrote %s', csv_path.name)
            draw_page(pdf, adata, version, rows, fc_df)
            logging.info('Page done: %s', version)

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — Fold-change curated dot plot (3 cell types)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "New version of the curated dot plot:\n"
          "  • CONSTANT dot size (fraction expressing NOT encoded).\n"
          "  • COLOUR = log2 fold change vs a reference group; diverging\n"
          f"    blue-white-red, centred 0, clipped ±{CLIP:g}.\n"
          "  • THREE cell-type panels: Excitatory / Inhibitory / Astrocyte.\n"
          "  • Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).\n\n"
          "Five columns per panel:\n"
          "  1. MCT1i        = H2O_MCT1i   vs H2O_veh\n"
          "  2. EtOH         = EtOH_veh    vs H2O_veh\n"
          "  3. EtOH+MCT1i   = EtOH_MCT1i  vs H2O_veh\n"
          "  4. Chronic EtOH = ChronicEtOH vs H2O_veh\n"
          "  5. EtOH+MCT1i vs EtOH = EtOH_MCT1i vs EtOH_veh  (MCT1i rescue: does\n"
          "     blocking monocarboxylate uptake reverse the ethanol effect? —\n"
          "     shown after a dashed divider, different reference). MAT2A/APP excluded.\n\n"
          "Fold change is a Slides-A+B pooled pseudobulk ratio per cell type:\n"
          "  log2FC = log2((mean_norm_expr[test]+1)/(mean_norm_expr[ref]+1)),\n"
          "mean_norm_expr = mean of linear (expm1 of log-norm) expression over all\n"
          "cells of that cell type in that group.\n\n"
          "Genes are grouped into 17 curated functional categories, followed by an\n"
          "'Additional user gene list' section (97 provided genes + GFP). Genes can\n"
          "repeat between the functional categories and the additional section.\n"
          "Blank cell = gene absent from panel or fewer than "
          f"{MIN_CELLS} cells in that group×cell-type;\n"
          "† after a gene name = not present in the Xenium panel.\n\n"
          "Caveat: Inhibitory neurons are sparse (especially in V2), so several\n"
          "Inhibitory cells may be blank/under-powered — interpret with care.\n"
          "Statistics caveat: cell-level pseudobulk means; treat direction and\n"
          "log2FC magnitude as primary, not significance.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=len(versions) - 1)
    out_pdf = args.outdir / f'Xenium_FCcuratedDotplot_3CT_V1V2{clip_tag}_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + len(versions))

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    logging.info('Done.')


if __name__ == '__main__':
    main()
