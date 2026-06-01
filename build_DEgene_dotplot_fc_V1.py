#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled, V1 only — FOLD-CHANGE dot plot of the DE
genes used in the cross-cell-type heat maps.

Fold-change counterpart of build_DEgene_dotplot_abundance_V1.py and a direct
dot-plot rendering of build_heatmaps_crossCT.py. The cross-CT heat maps anchor
on one cell type's DE gene list (EtOH_veh vs H2O_veh: top 50 up + 50 down) and
show those SAME genes' log2FC across every cell type, one column per contrast.
Here we reuse exactly that gene selection AND exactly the heat maps' log2FC
values (read from the same long-form DE CSV), but draw a DOT PLOT:

  * DOT SIZE  = % of cells positive for the transcript (counts > 0) in the
    TEST group of that contrast x cell type (computed from the cell data).
  * COLOUR ("heat") = log2 fold change (from the DE table — identical to the
    heat maps). Diverging blue-white-red (RdBu_r), clipped to +/-CLIP
    (default 0.6).
  * THREE cell-type panels: Excitatory / Inhibitory / Astrocyte (each gene read
    straight across cell types, exactly like the cross-CT heat maps).
  * V1 only (broad-ROI punches).

Columns = the FOUR contrasts the heat maps used:
    EtOH vs H2O   (EtOH_veh   vs H2O_veh)
    MCT1i vs H2O  (H2O_MCT1i  vs H2O_veh)
    EtOH+MCT1i vs H2O   (EtOH_MCT1i vs H2O_veh)
    EtOH+MCT1i vs EtOH  (EtOH_MCT1i vs EtOH_veh)

Gene (row) selection — identical to the heat maps' anchor logic:
    in the anchor cell type's EtOH_veh-vs-H2O_veh contrast,
    |log2FC| >= LFC_THRESH AND padj < PADJ_THRESH, then top TOP_N up + TOP_N down
    by EtOH log2FC. Anchors: Excitatory and Astrocyte (Inhibitory ~0 sig genes).

Cell typing follows the master build (argmax marker-module score; Excit/Inhib
split within neurons). All paths via CLI. Vector PDFs, editable text
(pdf.fonttype=42). One source CSV per anchor written alongside the figure.
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
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ----------------------------- parameters ----------------------------------
# (test, reference, column label) — must match the cross-CT heat maps' contrasts
CONTRASTS = [
    ('EtOH_veh',   'H2O_veh',  'EtOH\nvs H2O'),
    ('H2O_MCT1i',  'H2O_veh',  'MCT1i\nvs H2O'),
    ('EtOH_MCT1i', 'H2O_veh',  'EtOH+MCT1i\nvs H2O'),
    ('EtOH_MCT1i', 'EtOH_veh', 'E+M\nvs EtOH'),
]
GROUPS = sorted(set([t for t, _r, _l in CONTRASTS] + [r for _t, r, _l in CONTRASTS]))
PANELS = ['Excitatory', 'Inhibitory', 'Astrocyte']      # 'subtype' column
ANCHORS = ['Excitatory', 'Astrocyte']                   # one page each

# DE gene-selection (must match the cross-CT heat maps)
SEL_TEST, SEL_REF = 'EtOH_veh', 'H2O_veh'
LFC_THRESH = 0.5
PADJ_THRESH = 1e-3
TOP_N = 50

MIN_CELLS = 10          # min cells per (group x cell type) to plot a value
CMAP = 'RdBu_r'         # diverging "heat" for log2 fold change
DEFAULT_CLIP = 0.6      # log2FC colour clip (+/-)
SIZE_MIN, SIZE_MAX = 6.0, 190.0
LEGEND_FRACS = [0.25, 0.5, 0.75, 1.0]

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
    X = adata[:, gene].layers['counts']
    return np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()


def disp(g):
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


def size_of(frac):
    return SIZE_MIN + np.clip(frac, 0.0, 1.0) * (SIZE_MAX - SIZE_MIN)


# ----------------------------- DE gene selection ----------------------------
def anchor_genes(de, anchor):
    """Top TOP_N up + TOP_N down EtOH-vs-H2O genes in the anchor cell type
    (|log2FC| >= LFC_THRESH AND padj < PADJ_THRESH). Returns (up, down)."""
    sub = de[(de['subtype'] == anchor) & (de['test'] == SEL_TEST) &
             (de['reference'] == SEL_REF)]
    if sub.empty:
        return [], []
    lfc = sub.set_index('gene')['logfc']
    padj = sub.set_index('gene')['padj']
    sig = (lfc.abs() >= LFC_THRESH) & (padj < PADJ_THRESH)
    sel = lfc[sig]
    up = sel[sel > 0].sort_values(ascending=False).head(TOP_N)
    dn = sel[sel < 0].sort_values(ascending=True).head(TOP_N)
    return list(up.index), list(dn.index)


def anchor_rows(anchor, up, dn):
    """Ordered (section, gene): strongest up at top, strongest down at bottom."""
    rows = [(f'Up in EtOH ({anchor})', g) for g in up]
    rows += [(f'Down in EtOH ({anchor})', g) for g in dn[::-1]]
    return rows


# ----------------------------- % positive (dot size) ------------------------
def frac_table(adata, genes):
    """Long-form per (gene, cell type, group): frac_pos (dot size) and n_cells.
    Only the four groups used by the contrasts are computed."""
    sub = adata[(adata.obs['subtype'].astype(str).isin(PANELS)) &
                (adata.obs['group'].astype(str).isin(GROUPS))].copy()
    var = set(sub.var_names)
    present = [g for g in genes if g in var]

    grp = sub.obs['group'].astype(str).values
    ctv = sub.obs['subtype'].astype(str).values
    n_cells = (pd.DataFrame({'ct': ctv, 'group': grp})
               .groupby(['ct', 'group']).size())

    frac = {}
    for gene in present:
        c = counts_col(sub, gene)
        pos = (c > 0).astype(float)
        f = (pd.DataFrame({'pos': pos, 'ct': ctv, 'group': grp})
             .groupby(['ct', 'group'])['pos'].mean())
        for key in f.index:
            frac[(gene, *key)] = f[key]

    out = []
    for gene in genes:
        in_panel = gene in var
        for ct in PANELS:
            for gkey in GROUPS:
                n = int(n_cells.get((ct, gkey), 0))
                fv = frac.get((gene, ct, gkey), np.nan) if in_panel else np.nan
                out.append({'gene': gene, 'in_panel': in_panel, 'cell_type': ct,
                            'group': gkey, 'n_cells': n,
                            'frac_pos': (float(fv) if not pd.isna(fv) else np.nan),
                            'pct_pos': (float(fv) * 100.0 if not pd.isna(fv) else np.nan)})
    return pd.DataFrame(out)


def logfc_lookup(de):
    """(gene, subtype, test, reference) -> log2FC, from the DE table the heat
    maps used (so colours are identical to the heat maps)."""
    keep = de[de['test'].isin([t for t, _r, _l in CONTRASTS]) &
              de['reference'].isin([r for _t, r, _l in CONTRASTS])]
    lut = {}
    for r in keep.itertuples(index=False):
        lut[(r.gene, r.subtype, r.test, r.reference)] = float(r.logfc)
    return lut


# ----------------------------- drawing --------------------------------------
def build_panel_df(genes, frac_df, lfc_lut):
    """Per (gene, cell_type, contrast): log2FC (colour), frac_pos (size),
    n_test_cells, n_ref_cells. Under-powered / missing blanked (NaN)."""
    fr = {(r.gene, r.cell_type, r.group): (r.frac_pos, r.n_cells, r.in_panel)
          for r in frac_df.itertuples(index=False)}
    rows = []
    for gene in genes:
        for ct in PANELS:
            for test, ref, lab in CONTRASTS:
                fv, n_test, in_panel = fr.get((gene, ct, test), (np.nan, 0, False))
                _rf, n_ref, _ip = fr.get((gene, ct, ref), (np.nan, 0, False))
                lfc = lfc_lut.get((gene, ct, test, ref), np.nan)
                ok = in_panel and n_test >= MIN_CELLS and n_ref >= MIN_CELLS
                rows.append({
                    'gene': gene, 'in_panel': in_panel, 'cell_type': ct,
                    'contrast': lab.replace('\n', ' '), 'test': test, 'reference': ref,
                    'n_test_cells': int(n_test), 'n_ref_cells': int(n_ref),
                    'frac_pos_test': (float(fv) if ok and not pd.isna(fv) else np.nan),
                    'pct_pos_test': (float(fv) * 100.0 if ok and not pd.isna(fv) else np.nan),
                    'log2fc': (float(lfc) if ok and not pd.isna(lfc) else np.nan),
                })
    return pd.DataFrame(rows)


def draw_page(pdf, anchor, rows, df, panel_counts, norm, cmap, clip):
    n_rows = len(rows)
    if n_rows == 0:
        return
    y_of = {(sec, g): i for i, (sec, g) in enumerate(rows)}

    sec_order, sec_start = [], {}
    for i, (sec, _g) in enumerate(rows):
        if sec not in sec_start:
            sec_start[sec] = i
            sec_order.append(sec)
    sec_bounds = [sec_start[s] for s in sec_order] + [n_rows]

    ncol = len(CONTRASTS)
    col_labels = [lab for _t, _r, lab in CONTRASTS]
    fig_w = 3.6 + len(PANELS) * (0.66 * ncol + 0.8)
    fig_h = max(7.0, 0.16 * n_rows + 2.4)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(fig_w, fig_h), sharey=True,
                             gridspec_kw={'wspace': 0.12})
    if len(PANELS) == 1:
        axes = [axes]

    for ax, ct in zip(axes, PANELS):
        d = df[df['cell_type'] == ct]
        lut = {(r.gene, r.contrast): (r.frac_pos_test, r.log2fc)
               for r in d.itertuples(index=False)}
        for (sec, gene), yi in y_of.items():
            for j, (_t, _r, lab) in enumerate(CONTRASTS):
                fv, lfc = lut.get((gene, lab.replace('\n', ' ')), (np.nan, np.nan))
                if pd.isna(fv) or pd.isna(lfc) or fv <= 0:
                    continue
                ax.scatter(j, yi, s=size_of(fv), c=[lfc], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.3, zorder=3)
        ax.set_xlim(-0.6, ncol - 0.4)
        ax.set_ylim(n_rows - 0.5, -0.5)
        ax.set_xticks(range(ncol))
        ax.set_xticklabels(col_labels, rotation=0, ha='center', va='top', fontsize=5.6)
        secax = ax.secondary_xaxis('top')
        secax.set_xticks(range(ncol))
        secax.set_xticklabels(col_labels, rotation=0, ha='center', va='bottom', fontsize=5.6)
        secax.tick_params(length=0)
        for sp in secax.spines.values():
            sp.set_visible(False)
        ax.set_title(f'{ct}\n(n={panel_counts.get(ct, 0):,} cells)',
                     fontsize=9, fontweight='bold', pad=24)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        for b in sec_bounds[1:-1]:
            ax.axhline(b - 0.5, color='0.4', lw=0.8, ls='--', zorder=1)
        ax.set_axisbelow(True)
        ax.grid(axis='x', color='0.92', lw=0.4, zorder=0)

    present_genes = set(df.loc[df['in_panel'], 'gene'])
    ylabels, ycolors = [], []
    for _s, g in rows:
        if g in present_genes:
            ylabels.append(disp(g)); ycolors.append('black')
        else:
            ylabels.append(disp(g) + ' †'); ycolors.append('0.6')
    fs = max(3.0, min(6.5, 560 / max(n_rows, 1)))
    axes[0].set_yticks(range(n_rows))
    axes[0].set_yticklabels(ylabels, fontsize=fs)
    for tick, col in zip(axes[0].get_yticklabels(), ycolors):
        tick.set_color(col)

    rax = axes[-1]
    for s, e, sec in zip(sec_bounds[:-1], sec_bounds[1:], sec_order):
        yc = (s + e - 1) / 2.0
        rax.text(ncol - 0.30, yc, sec, fontsize=7.0, va='center', ha='left',
                 rotation=90, clip_on=False, color='0.15', fontweight='bold')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.16, extend='both')
    cbar.set_label(f'log2 fold change (clipped ±{clip:g})', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    handles = [ax.scatter([], [], s=size_of(f), c='0.55', edgecolor='black',
                          linewidth=0.3) for f in LEGEND_FRACS]
    labels = [f'{int(f * 100)}%' for f in LEGEND_FRACS]
    fig.legend(handles, labels, title='% cells positive\n(test group; dot size)',
               loc='center left', bbox_to_anchor=(0.995, 0.30), frameon=False,
               labelspacing=1.5, handletextpad=0.9, borderpad=0.8,
               fontsize=7, title_fontsize=7)

    n_up = sum(1 for s, _ in rows if s.startswith('Up'))
    n_dn = n_rows - n_up
    n_absent = sum(1 for _s, g in rows if g not in present_genes)
    fig.suptitle(
        f'V1 — fold-change dot plot of {anchor}-anchored DE genes (cross-cell-type)\n'
        f'rows = top {TOP_N} up + {TOP_N} down by EtOH log2FC in {anchor} '
        f'(|log2FC|≥{LFC_THRESH:g}, padj<{PADJ_THRESH:g}); shown n={n_up} up / {n_dn} down\n'
        f'dot SIZE = % cells positive (test group); COLOUR = log2FC (clipped ±{clip:g}); '
        f'3 cell types; columns = 4 contrasts; † not in panel (n={n_absent})',
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
    ap.add_argument('--de-csv', required=True, type=Path,
                    help='V1 long-form DE table the heat maps used (gene, logfc, '
                         'padj, subtype, test, reference, ...).')
    ap.add_argument('--clip', type=float, default=DEFAULT_CLIP,
                    help=f'log2FC colour clip (+/-); default {DEFAULT_CLIP:g}.')
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)
    clip = float(args.clip)
    clip_tag = '' if abs(clip - DEFAULT_CLIP) < 1e-9 else f"_clip{('%g' % clip).replace('.', 'p')}"

    de = pd.read_csv(args.de_csv)
    logging.info('DE table: %s rows; subtypes=%s', len(de), sorted(de['subtype'].unique()))

    anchor_data = {}
    all_genes = []
    for anchor in ANCHORS:
        up, dn = anchor_genes(de, anchor)
        rows = anchor_rows(anchor, up, dn)
        anchor_data[anchor] = rows
        all_genes.extend([g for _s, g in rows])
        logging.info('Anchor %s: %d up + %d down = %d genes', anchor, len(up), len(dn), len(rows))
    all_genes = list(dict.fromkeys(all_genes))   # unique, order-preserving

    lfc_lut = logfc_lookup(de)

    logging.info('=== Processing V1 (pool A+B) ===')
    adata = process_version(args.slide_a_dir, args.v1_slide_a_ann,
                            args.slide_b_dir, args.v1_slide_b_ann)
    panel_counts = (adata.obs[adata.obs['group'].astype(str).isin(GROUPS)]
                    ['subtype'].astype(str).value_counts().to_dict())
    logging.info('V1 subtype composition (contrast groups): %s', panel_counts)

    frac_df = frac_table(adata, all_genes)
    panel_df_all = build_panel_df(all_genes, frac_df, lfc_lut)

    norm = Normalize(vmin=-clip, vmax=clip, clip=True)
    cmap = plt.get_cmap(CMAP)
    logging.info('Colour scale: log2FC clipped ±%.3g (RdBu_r)', clip)

    tmp = args.outdir / '_tmp_de_fc_dot'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'
    with PdfPages(pages_pdf) as pdf:
        for anchor in ANCHORS:
            rows = anchor_data[anchor]
            genes = [g for _s, g in rows]
            df = panel_df_all[panel_df_all['gene'].isin(genes)].copy()
            sec_of = {g: s for s, g in rows}
            df['section'] = df['gene'].map(sec_of)
            csv_path = (args.outdir /
                        f'SlidesAB_DEgeneDotplot_FC_V1_anchor{anchor}{clip_tag}_{today}.csv')
            df.to_csv(csv_path, index=False)
            logging.info('Wrote %s', csv_path.name)
            draw_page(pdf, anchor, rows, df, panel_counts, norm, cmap, clip)
            logging.info('Page done: anchor=%s', anchor)

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — fold-change dot plot of cross-CT heat-map DE genes (V1)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "Dot-plot rendering of the cross-cell-type heat maps. The SAME anchor\n"
          "gene selection AND the SAME log2 fold changes (from the DE table) are\n"
          "reused; each dot encodes:\n"
          "  • DOT SIZE  = % of cells positive (counts > 0) in the TEST group of\n"
          "    the contrast.\n"
          f"  • COLOUR ('heat') = log2 fold change, diverging blue-white-red,\n"
          f"    clipped to ±{clip:g}.\n"
          "  • THREE cell-type panels: Excitatory / Inhibitory / Astrocyte (each\n"
          "    gene read straight across cell types, as in the heat maps).\n"
          "  • V1 only (broad-ROI punches).\n\n"
          "Columns = the four CONTRASTS used by the heat maps:\n"
          "  EtOH vs H2O · MCT1i vs H2O · EtOH+MCT1i vs H2O · EtOH+MCT1i vs EtOH.\n\n"
          "Gene (row) selection — identical to the heat maps' anchor logic:\n"
          f"  in the anchor cell type's {SEL_TEST} vs {SEL_REF} contrast,\n"
          f"  |log2FC| >= {LFC_THRESH:g} AND padj < {PADJ_THRESH:g}; then top {TOP_N} up + {TOP_N} down\n"
          "  by EtOH log2FC. One page per anchor: Excitatory, Astrocyte\n"
          "  (Inhibitory omitted — ~0 significant genes).\n"
          "Dashed line separates up- from down-regulated genes.\n"
          "† after a gene = not present in the Xenium panel after filtering;\n"
          f"blank cell = fewer than {MIN_CELLS} cells in test or reference group.\n\n"
          "Statistics caveat: cell-level Wilcoxon DE — adj p-values inflated by\n"
          "pseudoreplication; treat direction and log2FC magnitude as primary.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=len(ANCHORS) - 1)
    out_pdf = args.outdir / f'Xenium_DEgeneDotplot_FC_V1{clip_tag}_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + len(ANCHORS))

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    logging.info('Done.')


if __name__ == '__main__':
    main()
