#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled — classic dot plot (SIZE + HEAT), 3 cell types.

Third companion to build_fc_curated_dotplot_3CT_V1V2.py (fold change) and
build_pctpos_curated_dotplot_3CT_V1V2.py (% positive). Same curated gene layout
and 3 cell-type panels, but the dot encodes TWO quantities at once — the classic
single-cell "dot plot":

  * DOT SIZE  = % of cells positive for that transcript (fraction with detected
    counts > 0) in that group x cell type. Bigger dot = more cells express it.
  * COLOUR ("heat") = mean number of TRANSCRIPTS per cell (raw counts averaged
    over ALL cells of that group x cell type, zeros included). Log colour scale,
    SHARED across V1 and V2 so the two pages are directly comparable.
  * THREE cell-type panels: Excitatory / Inhibitory / Astrocyte.
  * Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).

Both quantities are per-group (not contrasts), so columns are the five alcohol
GROUPS themselves:
    H2O_veh, H2O_MCT1i, EtOH_veh, EtOH_MCT1i, ChronicEtOH
MAT2A / APP groups are excluded.

A group x cell-type cell-count guard (>= MIN_CELLS) blanks under-powered cells
(notably Inhibitory in V2).

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
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ----------------------------- parameters ----------------------------------
# Columns (left -> right): (group key, two-line label). One column per group.
GROUP_COLS = [
    ('H2O_veh',     'H2O\nveh'),
    ('H2O_MCT1i',   'H2O\nMCT1i'),
    ('EtOH_veh',    'EtOH\nveh'),
    ('EtOH_MCT1i',  'EtOH\nMCT1i'),
    ('ChronicEtOH', 'Chronic\nEtOH'),
]
GROUPS = [g for g, _ in GROUP_COLS]
PANELS = ['Excitatory', 'Inhibitory', 'Astrocyte']     # 'subtype' column
MIN_CELLS = 10          # min cells per (group x cell type) to plot a value
CMAP = 'viridis'        # sequential "heat" for mean transcript count
SIZE_MIN, SIZE_MAX = 6.0, 190.0   # marker-area range for 0% .. 100% positive
LEGEND_FRACS = [0.25, 0.5, 0.75, 1.0]   # size-legend reference detection rates

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
    adata.layers['counts'] = adata.X.copy()     # raw transcript counts (pre-norm)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return cell_type_data(adata)


def counts_col(adata, gene):
    """Raw transcript counts for one gene (from the stashed 'counts' layer)."""
    X = adata[:, gene].layers['counts']
    return np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()


def disp(g):
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


def size_of(frac):
    """Map detection fraction (0..1) to marker area."""
    return SIZE_MIN + np.clip(frac, 0.0, 1.0) * (SIZE_MAX - SIZE_MIN)


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


# ----------------------------- size + heat table ----------------------------
def size_heat_table(adata, rows):
    """Return long-form DataFrame per (gene, cell type, group) with:
        frac_pos        -> dot SIZE  (fraction of cells with counts > 0)
        mean_count      -> dot HEAT  (mean raw transcripts/cell, all cells)
        mean_count_pos  -> reference (mean raw transcripts among positive cells)
    plus n_cells. Under-powered / absent cells are blanked (NaN)."""
    sub = adata[(adata.obs['subtype'].astype(str).isin(PANELS)) &
                (adata.obs['group'].astype(str).isin(GROUPS))].copy()
    var = set(sub.var_names)
    genes_needed = sorted({g for _, g in rows})
    present = [g for g in genes_needed if g in var]

    grp = sub.obs['group'].astype(str).values
    ctv = sub.obs['subtype'].astype(str).values
    n_cells = (pd.DataFrame({'ct': ctv, 'group': grp})
               .groupby(['ct', 'group']).size())

    frac, mean_all, mean_pos = {}, {}, {}
    for gene in present:
        c = counts_col(sub, gene)                       # raw counts
        pos = (c > 0).astype(float)
        d = pd.DataFrame({'c': c, 'pos': pos, 'ct': ctv, 'group': grp})
        g = d.groupby(['ct', 'group'])
        f = g['pos'].mean()
        ma = g['c'].mean()
        # mean among positive cells = sum(c) / n_pos
        sums = g['c'].sum()
        npos = g['pos'].sum()
        for key in f.index:
            frac[(gene, *key)] = f[key]
            mean_all[(gene, *key)] = ma[key]
            mean_pos[(gene, *key)] = (sums[key] / npos[key]) if npos[key] > 0 else np.nan

    out = []
    for section, gene in rows:
        in_panel = gene in var
        for ct in PANELS:
            for gkey, _lab in GROUP_COLS:
                n = int(n_cells.get((ct, gkey), 0))
                if (not in_panel) or n < MIN_CELLS:
                    fv = mv = mp = np.nan
                else:
                    fv = frac.get((gene, ct, gkey), np.nan)
                    mv = mean_all.get((gene, ct, gkey), np.nan)
                    mp = mean_pos.get((gene, ct, gkey), np.nan)
                out.append({'section': section, 'gene': gene, 'in_panel': in_panel,
                            'cell_type': ct, 'group': gkey, 'n_cells': n,
                            'frac_pos': (float(fv) if not pd.isna(fv) else np.nan),
                            'pct_pos': (float(fv) * 100.0 if not pd.isna(fv) else np.nan),
                            'mean_count': (float(mv) if not pd.isna(mv) else np.nan),
                            'mean_count_pos': (float(mp) if not pd.isna(mp) else np.nan)})
    return pd.DataFrame(out)


# ----------------------------- drawing --------------------------------------
def draw_page(pdf, version, rows, df, panel_counts, norm, cmap):
    n_rows = len(rows)
    y_of = {(sec, g): i for i, (sec, g) in enumerate(rows)}  # 0 = top

    sec_order, sec_start = [], {}
    for i, (sec, _g) in enumerate(rows):
        if sec not in sec_start:
            sec_start[sec] = i
            sec_order.append(sec)
    sec_bounds = [sec_start[s] for s in sec_order] + [n_rows]

    ncol = len(GROUP_COLS)

    fig_w = 3.6 + len(PANELS) * (0.62 * ncol + 0.8)
    fig_h = max(8.0, 0.135 * n_rows + 2.4)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(fig_w, fig_h), sharey=True,
                             gridspec_kw={'wspace': 0.12})
    if len(PANELS) == 1:
        axes = [axes]

    for ax, ct in zip(axes, PANELS):
        d = df[df['cell_type'] == ct]
        lut = {(r.gene, r.section, r.group): (r.frac_pos, r.mean_count)
               for r in d.itertuples(index=False)}
        for (sec, gene), yi in y_of.items():
            for j, (gkey, _lab) in enumerate(GROUP_COLS):
                fv, mv = lut.get((gene, sec, gkey), (np.nan, np.nan))
                if pd.isna(fv) or pd.isna(mv) or fv <= 0:
                    continue
                ax.scatter(j, yi, s=size_of(fv), c=[mv], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.3, zorder=3)
        ax.set_xlim(-0.6, ncol - 0.4)
        ax.set_ylim(n_rows - 0.5, -0.5)   # row 0 at top
        col_labels = [lab for _g, lab in GROUP_COLS]
        ax.set_xticks(range(ncol))
        ax.set_xticklabels(col_labels, rotation=0, ha='center', va='top', fontsize=6.0)
        secax = ax.secondary_xaxis('top')
        secax.set_xticks(range(ncol))
        secax.set_xticklabels(col_labels, rotation=0, ha='center', va='bottom', fontsize=6.0)
        secax.tick_params(length=0)
        for sp in secax.spines.values():
            sp.set_visible(False)
        ax.set_title(f'{ct}\n(n={panel_counts.get(ct, 0):,} cells)',
                     fontsize=9, fontweight='bold', pad=24)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        for b in sec_bounds[1:-1]:
            ax.axhline(b - 0.5, color='0.6', lw=0.6, zorder=1)
        ax.set_axisbelow(True)
        ax.grid(axis='x', color='0.92', lw=0.4, zorder=0)

    # y labels (gene names) on leftmost axis; mark genes absent from panel
    present_genes = set(df.loc[df['in_panel'], 'gene'])
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

    # heat colorbar (mean transcripts/cell, log scale)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.16)
    cbar.set_label('mean transcripts per cell (log scale)', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # size legend (% cells positive) — proxy handles
    handles = [ax.scatter([], [], s=size_of(f), c='0.55', edgecolor='black',
                          linewidth=0.3) for f in LEGEND_FRACS]
    labels = [f'{int(f * 100)}%' for f in LEGEND_FRACS]
    fig.legend(handles, labels, title='% cells positive\n(dot size)',
               loc='center left', bbox_to_anchor=(0.995, 0.30), frameon=False,
               labelspacing=1.5, handletextpad=0.9, borderpad=0.8,
               fontsize=7, title_fontsize=7)

    n_absent = sum(1 for _s, g in rows if g not in present_genes)
    fig.suptitle(
        f'Slides A+B pooled — dot plot (size = % positive, heat = transcripts) — {version}\n'
        f'dot SIZE = % of cells expressing the transcript; '
        f'COLOUR = mean transcripts/cell (log, all cells); 3 cell types; columns = 5 alcohol groups\n'
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
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = build_rows()
    logging.info('Total rows: %d (%d functional + %d extra-section)',
                 len(rows), sum(1 for s, _ in rows if s != EXTRA_SECTION),
                 sum(1 for s, _ in rows if s == EXTRA_SECTION))

    tmp = args.outdir / '_tmp_countheat_dot'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    versions = {
        'V1': (args.v1_slide_a_ann, args.v1_slide_b_ann),
        'V2': (args.v2_slide_a_ann, args.v2_slide_b_ann),
    }

    # Pass 1: build per-version tables + panel counts; collect heat values for a
    # SHARED log color scale across V1 and V2.
    results = []
    all_means = []
    for version, (a_ann, b_ann) in versions.items():
        logging.info('=== Processing %s (pool A+B) ===', version)
        adata = process_version(args.slide_a_dir, a_ann, args.slide_b_dir, b_ann)
        panel_counts = (adata.obs[adata.obs['group'].astype(str).isin(GROUPS)]
                        ['subtype'].astype(str).value_counts().to_dict())
        logging.info('%s subtype composition (alcohol groups): %s', version, panel_counts)
        df = size_heat_table(adata, rows)
        csv_path = args.outdir / f'SlidesAB_CountHeatCuratedDotplot_3CT_{version}_{today}.csv'
        df.to_csv(csv_path, index=False)
        logging.info('Wrote %s', csv_path.name)
        results.append((version, df, panel_counts))
        all_means.extend(df.loc[df['mean_count'] > 0, 'mean_count'].tolist())
        del adata

    arr = np.asarray([v for v in all_means if v > 0], dtype=float)
    vmin = max(float(np.percentile(arr, 1.0)), 1e-2)
    vmax = float(np.percentile(arr, 99.5))
    if vmax <= vmin:
        vmax = vmin * 10
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)
    cmap = plt.get_cmap(CMAP)
    logging.info('Shared heat scale (mean transcripts/cell): vmin=%.3g vmax=%.3g', vmin, vmax)

    # Pass 2: draw pages with the shared color scale.
    with PdfPages(pages_pdf) as pdf:
        for version, df, panel_counts in results:
            draw_page(pdf, version, rows, df, panel_counts, norm, cmap)
            logging.info('Page done: %s', version)

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — dot plot: size = % positive, heat = transcripts (3 cell types)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "Third companion to the curated dot plots (fold change; % positive).\n"
          "Same curated gene layout and 3 cell-type panels, but each dot encodes\n"
          "TWO quantities at once — the classic single-cell dot plot:\n"
          "  • DOT SIZE  = % of cells positive for the transcript (counts > 0);\n"
          "    bigger dot = more cells express it (0-100%).\n"
          "  • COLOUR ('heat') = mean number of TRANSCRIPTS per cell (raw counts\n"
          "    averaged over ALL cells of that group x cell type, zeros included).\n"
          "    Log colour scale, SHARED across V1 and V2 so pages are comparable.\n"
          "  • THREE cell-type panels: Excitatory / Inhibitory / Astrocyte.\n"
          "  • Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).\n\n"
          "Columns = the five alcohol GROUPS themselves (both quantities are\n"
          "per-group, not contrasts):\n"
          "  H2O_veh, H2O_MCT1i, EtOH_veh, EtOH_MCT1i, ChronicEtOH.\n"
          "MAT2A / APP groups excluded.\n\n"
          "Genes are grouped into 17 curated functional categories, followed by an\n"
          "'Additional user gene list' section (97 provided genes + GFP). Genes can\n"
          "repeat between the functional categories and the additional section.\n"
          "Blank cell = gene absent from panel or fewer than "
          f"{MIN_CELLS} cells in that group×cell-type;\n"
          "† after a gene name = not present in the Xenium panel.\n\n"
          "The source CSV also reports mean_count_pos (mean transcripts among\n"
          "positive cells only) for reference. Caveat: Inhibitory neurons are\n"
          "sparse (especially V2); detection/abundance also reflect probe\n"
          "efficiency — compare the same gene across groups, not different genes.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=len(versions) - 1)
    out_pdf = args.outdir / f'Xenium_CountHeatCuratedDotplot_3CT_V1V2_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + len(versions))

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    logging.info('Done.')


if __name__ == '__main__':
    main()
