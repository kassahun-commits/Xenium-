#!/usr/bin/env python3
"""
Xenium May 2026 — Neuron vs Astrocyte fold-change DOT PLOT, one figure per
SEGMENTATION ITERATION (whole-cell vs nucleus-only).

This generalises the alcohol-responsive dot plot
(make_alcohol_responsive_dotplot.py) so the deck can compare the SAME curated
gene panel across the two segmentation iterations on identical axes.

Encoding (per gene x cell type x contrast):
  * DOT SIZE  = % of cells positive (raw counts > 0) in the TEST group of the
                contrast (computed from each dataset's own counts).
  * COLOUR    = Wilcoxon log2 fold change (test / reference), read from the
                pre-computed marker-typed DE tables (NOT recomputed here):
                  - whole-cell : PI_Master_NvsA_V1_DE  (col 'logfc')
                  - nucleus    : Nucleus_WilcoxonDE    (col 'log2FoldChange')
                diverging blue-white-red (RdBu_r), clipped to [vmin, vmax] (±1).
  * PANELS    = Neuron, Astrocyte.
  * COLUMNS   = four contrasts, all vs H2O_veh (the only references the nucleus
                Wilcoxon table carries), so whole-cell and nucleus share axes:
                  EtOH_veh    vs H2O_veh   -> 'EtOH/con'
                  ChronicEtOH vs H2O_veh   -> 'chronic/con'
                  EtOH_MCT1i  vs H2O_veh   -> 'MCT1i+EtOH/con'
                  H2O_MCT1i   vs H2O_veh   -> 'MCT1i/con'

NOTE on FC convention: the original alcohol dot plot coloured by a
mean-expression ratio; here BOTH iterations are coloured by the Wilcoxon
log2FC so the deck's Wilcoxon thread is carried through and the two iterations
are directly comparable. Dot size (% positive) is unchanged.

MEWS rules: no hardcoded paths (CLI), editable vector text (pdf.fonttype=42),
source-data CSV alongside the figure. Outputs: PNG (300 dpi), PDF, CSV.
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
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# --------------------------- parameters -------------------------------------
PANELS = ['Neuron', 'Astrocyte']

# (test, reference, column label). All vs H2O_veh so whole-cell & nucleus match.
CONTRASTS = [
    ('EtOH_veh',    'H2O_veh', 'EtOH/con'),
    ('ChronicEtOH', 'H2O_veh', 'chronic/con'),
    ('EtOH_MCT1i',  'H2O_veh', 'MCT1i+EtOH/con'),
    ('H2O_MCT1i',   'H2O_veh', 'MCT1i/con'),
]
GROUPS = sorted(set([t for t, _r, _l in CONTRASTS] + [r for _t, r, _l in CONTRASTS]))

# curated gene panel, grouped top->bottom; group boundaries drive separators
GENE_GROUPS = [
    ('Gene of interest',   ['Grin2d']),
    ('Neuronal markers',   ['Rbfox3', 'Snap25', 'Syn1', 'Stmn2', 'Map2', 'Tubb3']),
    ('Astrocyte markers',  ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1']),
    ('Activity / receptors', ['Oprm1', 'Bdnf', 'Arc', 'Egr1', 'Fos']),
    ('Acetyl-CoA / chromatin', ['Acss1', 'Acss2', 'Acly', 'Kat2a', 'Kat2b', 'Ep300',
                                'Kat5', 'Crebbp', 'Hdac1', 'Hdac2', 'Hdac3']),
    ('Alcohol-related',    ['Gabra1', 'Gabra2', 'Gabbr1', 'Grin1', 'Grin2a', 'Grin2b',
                            'Gria1', 'Drd1', 'Drd2', 'Crhr1', 'Aldh2', 'Cyp2e1']),
]
DEFAULT_GENES = [g for _lab, gl in GENE_GROUPS for g in gl]

SIZE_MIN, SIZE_MAX = 12.0, 300.0
LEGEND_FRACS = [0.25, 0.5, 0.75, 1.0]


def size_of(frac):
    return SIZE_MIN + np.clip(frac, 0.0, 1.0) * (SIZE_MAX - SIZE_MIN)


def dense_col(M):
    return np.asarray(M.todense()).ravel() if hasattr(M, 'todense') else np.asarray(M).ravel()


# --------------------------- % positive from an AnnData ---------------------
def frac_positive(adata, genes, celltype_col='celltype'):
    """Per (gene, cell_type in PANELS, group in GROUPS): % cells with counts>0,
    mean (normalized/log) expression, and n cells per (cell_type, group).
    Uses layers['counts'] for detection; uses .X (log-normalized) for the mean
    expression that colours the baseline column."""
    keep = (adata.obs[celltype_col].astype(str).isin(PANELS) &
            adata.obs['group'].astype(str).isin(GROUPS))
    sub = adata[keep.values].copy()
    var = set(map(str, sub.var_names))
    present = [g for g in genes if g in var]

    ctv = sub.obs[celltype_col].astype(str).values
    grp = sub.obs['group'].astype(str).values
    keyidx = {(ct, g): np.where((ctv == ct) & (grp == g))[0]
              for ct in PANELS for g in GROUPS}
    n_cells = {k: int(v.size) for k, v in keyidx.items()}

    use_counts = 'counts' in sub.layers
    frac_pos, mean_expr = {}, {}
    for gene in present:
        c = dense_col(sub[:, gene].layers['counts'] if use_counts else sub[:, gene].X)
        x = dense_col(sub[:, gene].X)
        pos = (c > 0).astype(float)
        for k, idx in keyidx.items():
            if idx.size:
                frac_pos[(gene, *k)] = float(pos[idx].mean())
                mean_expr[(gene, *k)] = float(x[idx].mean())
    return frac_pos, mean_expr, n_cells, set(present)


# --------------------------- log2FC lookups ---------------------------------
def lfc_lookup_wholecell(csv):
    """PI_Master_NvsA_V1_DE: cols gene, logfc, subtype, test, reference."""
    d = pd.read_csv(csv)
    d = d[d['reference'] == 'H2O_veh']
    out = {}
    for r in d.itertuples(index=False):
        out[(str(r.gene), str(r.subtype), str(r.test))] = float(r.logfc)
    return out


def lfc_lookup_nucleus(csv):
    """Nucleus_WilcoxonDE: cols celltype, comparison (TEST_vs_H2O_veh),
    gene, log2FoldChange."""
    d = pd.read_csv(csv)
    out = {}
    for r in d.itertuples(index=False):
        comp = str(r.comparison)
        if not comp.endswith('_vs_H2O_veh'):
            continue
        test = comp[:-len('_vs_H2O_veh')]
        out[(str(r.gene), str(r.celltype), test)] = float(r.log2FoldChange)
    return out


# --------------------------- build table ------------------------------------
BASELINE_GROUP = 'H2O_veh'
BASELINE_LABEL = 'H2O_veh\nbaseline'


def build_table(genes, frac_pos, mean_expr, n_cells, present, lfc_lut,
                min_cells, segment=None, contrasts=CONTRASTS):
    rows = []
    for gene in genes:
        in_panel = gene in present
        for ct in PANELS:
            # --- baseline column: H2O_veh control (size=%pos, colour=mean expr) ---
            n_b = n_cells.get((ct, BASELINE_GROUP), 0)
            fb = frac_pos.get((gene, ct, BASELINE_GROUP), np.nan)
            mb = mean_expr.get((gene, ct, BASELINE_GROUP), np.nan)
            okb = in_panel and n_b >= min_cells
            rows.append({
                'segment': segment, 'gene': gene, 'in_panel': in_panel,
                'cell_type': ct, 'contrast': BASELINE_LABEL, 'test': BASELINE_GROUP,
                'reference': '', 'n_test_cells': int(n_b),
                'frac_pos_test': (float(fb) if okb and not pd.isna(fb) else np.nan),
                'pct_pos_test': (float(fb) * 100.0 if okb and not pd.isna(fb) else np.nan),
                'log2FC': np.nan,
                'mean_expr': (float(mb) if okb and not pd.isna(mb) else np.nan),
            })
            # --- fold-change columns (colour = Wilcoxon log2FC) ---
            for test, ref, lab in contrasts:
                n_t = n_cells.get((ct, test), 0)
                fv = frac_pos.get((gene, ct, test), np.nan)
                lfc = lfc_lut.get((gene, ct, test), np.nan)
                ok = in_panel and n_t >= min_cells
                rows.append({
                    'segment': segment, 'gene': gene, 'in_panel': in_panel,
                    'cell_type': ct, 'contrast': lab, 'test': test,
                    'reference': ref, 'n_test_cells': int(n_t),
                    'frac_pos_test': (float(fv) if ok and not pd.isna(fv) else np.nan),
                    'pct_pos_test': (float(fv) * 100.0 if ok and not pd.isna(fv) else np.nan),
                    'log2FC': (float(lfc) if ok and not pd.isna(lfc) else np.nan),
                    'mean_expr': np.nan,
                })
    return pd.DataFrame(rows)


# --------------------------- drawing ----------------------------------------
# baseline column geometry: the H2O_veh control sits to the LEFT of the
# fold-change columns, separated by a small gap + divider line.
BASE_X = -1.5          # x position of the baseline column
BASE_GAP = -0.75       # x of the divider between baseline and FC columns


def draw(df, genes, panel_counts, norm, cmap, vmin, vmax, present, title,
         out_png, out_pdf, base_norm, base_cmap, gene_fs=8.5, col_fs=8.0,
         panel_fs=11.0, size_scale=0.8):
    n = len(genes)
    y_of = {g: i for i, g in enumerate(genes)}
    ncol = len(CONTRASTS)
    col_labels = [lab for _t, _r, lab in CONTRASTS]

    def sz(frac):
        return size_of(frac) * size_scale

    fig_w = 3.2 + len(PANELS) * (0.52 * (ncol + 1) + 0.6)
    fig_h = max(6.0, 0.32 * n + 2.4)
    fig, axes = plt.subplots(1, len(PANELS), figsize=(fig_w, fig_h), sharey=True,
                             gridspec_kw={'wspace': 0.16})
    if len(PANELS) == 1:
        axes = [axes]

    # group separator y positions (between consecutive gene-group blocks)
    sep_after = []
    acc = 0
    for _lab, gl in GENE_GROUPS:
        present_in_group = [g for g in gl if g in genes]
        acc += len(present_in_group)
        sep_after.append(acc - 0.5)
    sep_after = sep_after[:-1]  # no line after last block

    for ax, ct in zip(axes, PANELS):
        d = df[df['cell_type'] == ct]
        lut = {(r.gene, r.contrast): (r.frac_pos_test, r.log2FC)
               for r in d.itertuples(index=False)}
        blut = {r.gene: (r.frac_pos_test, r.mean_expr)
                for r in d[d['contrast'] == BASELINE_LABEL].itertuples(index=False)}
        for gene, yi in y_of.items():
            # baseline dot: colour = mean expression in H2O_veh
            fb, me = blut.get(gene, (np.nan, np.nan))
            if not (pd.isna(fb) or pd.isna(me) or fb <= 0):
                ax.scatter(BASE_X, yi, s=sz(fb), c=[me], cmap=base_cmap,
                           norm=base_norm, edgecolor='black', linewidth=0.4,
                           zorder=3)
            # fold-change dots
            for j, (_t, _r, lab) in enumerate(CONTRASTS):
                fv, lfc = lut.get((gene, lab), (np.nan, np.nan))
                if pd.isna(fv) or pd.isna(lfc) or fv <= 0:
                    continue
                ax.scatter(j, yi, s=sz(fv), c=[lfc], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.4, zorder=3)
        ax.axvline(BASE_GAP, color='0.7', lw=1.0, linestyle=(0, (4, 3)), zorder=1)
        ax.set_xlim(BASE_X - 0.6, ncol - 0.4)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks([BASE_X] + list(range(ncol)))
        ax.set_xticklabels([BASELINE_LABEL] + col_labels, rotation=45,
                           ha='right', va='top', fontsize=col_fs)
        ax.set_title(f'{ct}\n(n={panel_counts.get(ct, 0):,})',
                     fontsize=panel_fs, fontweight='bold', pad=10)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(color='0.93', lw=0.5, zorder=0)
        for ys in sep_after:
            ax.axhline(ys, color='0.80', lw=0.8, zorder=1)

    ylabels, ycolors = [], []
    for g in genes:
        if g in present:
            ylabels.append(g); ycolors.append('black')
        else:
            ylabels.append(g + ' †'); ycolors.append('0.6')
    axes[0].set_yticks(range(n))
    axes[0].set_yticklabels(ylabels, fontsize=gene_fs, fontstyle='italic')
    for tick, col in zip(axes[0].get_yticklabels(), ycolors):
        tick.set_color(col)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.024, pad=0.04, extend='both')
    cbar.set_label('Wilcoxon log2 fold change (test / H2O_veh)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    smb = plt.cm.ScalarMappable(cmap=base_cmap, norm=base_norm); smb.set_array([])
    cbar2 = fig.colorbar(smb, ax=axes, fraction=0.024, pad=0.10, extend='max')
    cbar2.set_label('baseline mean expr (H2O_veh, log-norm)', fontsize=9)
    cbar2.ax.tick_params(labelsize=8)

    handles = [ax.scatter([], [], s=sz(f), c='0.55', edgecolor='black',
                          linewidth=0.4) for f in LEGEND_FRACS]
    labels = [f'{int(f * 100)}%' for f in LEGEND_FRACS]
    fig.legend(handles, labels, title='% cells\npositive',
               loc='center left', bbox_to_anchor=(0.995, 0.30), frameon=False,
               labelspacing=1.6, handletextpad=1.0, borderpad=0.8,
               fontsize=8, title_fontsize=8)

    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.998)
    fig.text(0.5, 0.972,
             'dot size = % cells positive · FC cols colour = Wilcoxon log2FC '
             f'clipped [{vmin:g}, {vmax:g}] · baseline col colour = mean expr '
             '· † = gene absent from panel',
             ha='center', va='top', fontsize=7.5, color='0.35')

    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


# --------------------------- combined drawing -------------------------------
# Generic 4-panel side-by-side figure: genes listed ONCE on the left, then
# two HALVES of two cell-type panels each. Sharing the gene axis lets the dots
# be large. Used for both:
#   * segmentation compare : whole-cell  vs nucleus-only   (same QC, same typing)
#   * pipeline    compare  : V1 DE QC    vs 10X workshop QC (same marker typing)
# In every case dot colour = Wilcoxon log2FC and the only thing that changes
# across the two halves is the named factor.
COMBINED_PANELS = [
    ('WholeCell', 'Neuron'), ('WholeCell', 'Astrocyte'),
    ('NucleusOnly', 'Neuron'), ('NucleusOnly', 'Astrocyte'),
]
# half spec = (title, subtitle, colour, k0, k1) where k0..k1 are panel indices
HALF_SPECS_SEG = [
    ('WHOLE-CELL', 'marker-gene typing · Wilcoxon', '#0E7C86', 0, 1),
    ('NUCLEUS-ONLY', 'marker-gene typing · Wilcoxon', '#8E44AD', 2, 3),
]


def draw_combined(df, genes, panel_counts, norm, cmap, vmin, vmax, present,
                  out_png, out_pdf, base_norm, base_cmap,
                  panels=COMBINED_PANELS, contrasts=CONTRASTS,
                  half_specs=HALF_SPECS_SEG, caption=None,
                  gene_fs=9.0, col_fs=8.0, panel_fs=11.5, size_scale=1.28):
    n = len(genes)
    y_of = {g: i for i, g in enumerate(genes)}
    ncol = len(contrasts)
    col_labels = [lab for _t, _r, lab in contrasts]
    npan = len(panels)

    def sz(frac):
        return size_of(frac) * size_scale

    fig = plt.figure(figsize=(16.0, 8.4))
    gs = fig.add_gridspec(1, npan, wspace=0.14,
                          left=0.075, right=0.88, top=0.88, bottom=0.12)
    axes = [fig.add_subplot(gs[0, k]) for k in range(npan)]

    sep_after, acc = [], 0
    for _lab, gl in GENE_GROUPS:
        acc += len([g for g in gl if g in genes])
        sep_after.append(acc - 0.5)
    sep_after = sep_after[:-1]

    for k, (ax, (seg, ct)) in enumerate(zip(axes, panels)):
        d = df[(df['segment'] == seg) & (df['cell_type'] == ct)]
        lut = {(r.gene, r.contrast): (r.frac_pos_test, r.log2FC)
               for r in d.itertuples(index=False)}
        blut = {r.gene: (r.frac_pos_test, r.mean_expr)
                for r in d[d['contrast'] == BASELINE_LABEL].itertuples(index=False)}
        for gene, yi in y_of.items():
            # baseline dot: colour = mean expression in H2O_veh
            fb, me = blut.get(gene, (np.nan, np.nan))
            if not (pd.isna(fb) or pd.isna(me) or fb <= 0):
                ax.scatter(BASE_X, yi, s=sz(fb), c=[me], cmap=base_cmap,
                           norm=base_norm, edgecolor='black', linewidth=0.4,
                           zorder=3)
            for j, (_t, _r, lab) in enumerate(contrasts):
                fv, lfc = lut.get((gene, lab), (np.nan, np.nan))
                if pd.isna(fv) or pd.isna(lfc) or fv <= 0:
                    continue
                ax.scatter(j, yi, s=sz(fv), c=[lfc], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.4, zorder=3)
        ax.axvline(BASE_GAP, color='0.7', lw=0.9, linestyle=(0, (3, 2)), zorder=1)
        ax.set_xlim(BASE_X - 0.6, ncol - 0.4)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks([BASE_X] + list(range(ncol)))
        ax.set_xticklabels([BASELINE_LABEL] + col_labels, rotation=45,
                           ha='right', va='top', fontsize=col_fs)
        ax.set_title(f'{ct}\n(n={panel_counts.get((seg, ct), 0):,})',
                     fontsize=panel_fs, fontweight='bold', pad=8)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(color='0.93', lw=0.5, zorder=0)
        for ys in sep_after:
            ax.axhline(ys, color='0.80', lw=0.8, zorder=1)
        if k == 0:
            ylabels, ycolors = [], []
            for g in genes:
                ylabels.append(g if g in present else g + ' †')
                ycolors.append('black' if g in present else '0.6')
            ax.set_yticks(range(n))
            ax.set_yticklabels(ylabels, fontsize=gene_fs, fontstyle='italic')
            for tick, col in zip(ax.get_yticklabels(), ycolors):
                tick.set_color(col)
        else:
            ax.set_yticks([])

    # half super-headers spanning each pair of panels (figure coords)
    for title, sub, color, k0, k1 in half_specs:
        b0 = axes[k0].get_position(); b1 = axes[k1].get_position()
        xc = (b0.x0 + b1.x1) / 2.0
        fig.text(xc, 0.965, title, ha='center', va='center',
                 fontsize=15, fontweight='bold', color=color)
        fig.text(xc, 0.928, sub, ha='center', va='center',
                 fontsize=8.5, color='0.4')
        fig.lines.append(plt.Line2D([b0.x0, b1.x1], [0.915, 0.915],
                         transform=fig.transFigure, color=color, lw=2.0))

    # divider between the two halves
    mid = npan // 2
    xmid = (axes[mid - 1].get_position().x1 + axes[mid].get_position().x0) / 2.0
    fig.lines.append(plt.Line2D([xmid, xmid], [0.13, 0.90],
                     transform=fig.transFigure, color='0.75', lw=1.0,
                     linestyle=(0, (4, 3))))

    # ---- right-hand legends: two colorbars stacked + size legend below ----
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cax1 = fig.add_axes([0.915, 0.58, 0.013, 0.29])
    cbar = fig.colorbar(sm, cax=cax1, extend='both')
    cbar.set_label('Wilcoxon log2FC\n(test / H2O_veh)', fontsize=8.5)
    cbar.ax.tick_params(labelsize=8)

    smb = plt.cm.ScalarMappable(cmap=base_cmap, norm=base_norm); smb.set_array([])
    cax2 = fig.add_axes([0.915, 0.21, 0.013, 0.29])
    cbar2 = fig.colorbar(smb, cax=cax2, extend='max')
    cbar2.set_label('baseline mean expr\n(H2O_veh, log-norm)', fontsize=8.5)
    cbar2.ax.tick_params(labelsize=8)

    handles = [axes[0].scatter([], [], s=sz(f), c='0.55', edgecolor='black',
                               linewidth=0.4) for f in LEGEND_FRACS]
    labels = [f'{int(f * 100)}%' for f in LEGEND_FRACS]
    fig.legend(handles, labels, title='% cells\npositive',
               loc='upper left', bbox_to_anchor=(0.905, 0.17), frameon=False,
               labelspacing=1.4, handletextpad=1.0, borderpad=0.6,
               fontsize=8, title_fontsize=8)

    if caption is None:
        caption = ('dot size = % cells positive · FC cols colour = Wilcoxon '
                   f'log2FC clipped [{vmin:g}, {vmax:g}] · leftmost "H2O_veh '
                   'baseline" col colour = mean expr · † = gene absent')
    fig.text(0.5, 0.025, caption, ha='center', va='center',
             fontsize=8, color='0.35')

    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


# panels / halves / contrasts for the V1-DE-vs-10X pipeline comparison
PIPELINE_PANELS = [
    ('V1_DE', 'Neuron'), ('V1_DE', 'Astrocyte'),
    ('TenX', 'Neuron'), ('TenX', 'Astrocyte'),
]
HALF_SPECS_PIPE = [
    ('V1 DE pipeline', 'QC min_counts 10 · min_cells 5 · marker-typed',
     '#0E7C86', 0, 1),
    ('10X workshop pipeline',
     'QC min_counts 20 · max 3405 · min_cells 100 · marker-typed',
     '#B7791F', 2, 3),
]
PIPELINE_CONTRASTS = [
    ('EtOH_veh',    'H2O_veh', 'EtOH/con'),
    ('ChronicEtOH', 'H2O_veh', 'chronic/con'),
]


# --------------------------- data loaders -----------------------------------
def load_wholecell(args):
    """Reuse the tracked V1 dot-plot loaders + marker celltyping verbatim, so
    the whole-cell cell typing matches the rest of the project exactly."""
    sys.path.insert(0, str(args.loader_script_dir))
    import build_DEgene_dotplot_fc_V1 as base  # noqa: E402
    adata = base.process_version(args.slide_a_dir, args.slide_a_ann,
                                 args.slide_b_dir, args.slide_b_ann)
    return adata  # obs has 'celltype' (Neuron/Astrocyte/...) + 'group', layers['counts']


def load_nucleus(args):
    adata = ad.read_h5ad(args.nucleus_h5ad)
    return adata


def load_10x(args):
    """The 10X-workshop-QC'd whole-cell matrix (min_counts 20 / max 3405 /
    min_cells 100), MARKER-gene typed with the SAME function as the V1-DE
    half, so the only difference between the two halves is the QC. Reuses the
    counts layer for % positive and the lognorm layer for scoring + baseline."""
    sys.path.insert(0, str(args.loader_script_dir))
    import build_DEgene_dotplot_fc_V1 as base  # noqa: E402
    a = ad.read_h5ad(args.tenx_h5ad)
    if 'lognorm' in a.layers:
        a.X = a.layers['lognorm'].copy()
    a = base.cell_type_data(a)   # adds obs['celltype'] Neuron/Astrocyte/...
    return a


def lfc_lookup_compare(test_csvs, method='marker-typed'):
    """V1_celltyping_compare_Wilcoxon_<test>_vs_H2O_veh CSVs
    (cols celltype, gene, logfc, padj, method). Pick one typing method;
    key (gene, celltype, test) -> logfc."""
    out = {}
    for test, path in test_csvs.items():
        d = pd.read_csv(path)
        d = d[d['method'].astype(str) == method]
        for r in d.itertuples(index=False):
            out[(str(r.gene), str(r.celltype), str(test))] = float(r.logfc)
    return out


# --------------------------- main -------------------------------------------
def fmt(v):
    return ('%g' % v).replace('-', 'neg').replace('.', 'p')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', required=True,
                    choices=['wholecell', 'nucleus', 'combined',
                             'pipeline_compare'])
    ap.add_argument('--genes', nargs='+', default=DEFAULT_GENES)
    ap.add_argument('--vmin', type=float, default=-1.0)
    ap.add_argument('--vmax', type=float, default=1.0)
    ap.add_argument('--cmap', default='RdBu_r')
    ap.add_argument('--baseline-cmap', default='viridis',
                    help='sequential cmap for the H2O_veh baseline column '
                         '(mean expression); kept distinct from the diverging '
                         'fold-change cmap')
    ap.add_argument('--min-cells', type=int, default=10)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--date', default=date.today().isoformat())
    # whole-cell inputs
    ap.add_argument('--slide-a-dir', type=Path)
    ap.add_argument('--slide-b-dir', type=Path)
    ap.add_argument('--slide-a-ann', type=Path)
    ap.add_argument('--slide-b-ann', type=Path)
    ap.add_argument('--wholecell-de-csv', type=Path,
                    help='PI_Master_NvsA_V1_DE csv (cols gene,logfc,subtype,test,reference)')
    ap.add_argument('--loader-script-dir', type=Path,
                    default=Path(__file__).resolve().parent,
                    help='dir containing build_DEgene_dotplot_fc_V1.py '
                         '(supplies the raw-data loaders); defaults to this '
                         "script's own directory")
    # nucleus inputs
    ap.add_argument('--nucleus-h5ad', type=Path)
    ap.add_argument('--nucleus-de-csv', type=Path,
                    help='Nucleus_WilcoxonDE csv (cols celltype,comparison,gene,log2FoldChange)')
    # 10X-pipeline inputs (mode=pipeline_compare)
    ap.add_argument('--tenx-h5ad', type=Path,
                    help='WholeCell_V1_10x_processed h5ad (10X workshop QC; '
                         'layers counts+lognorm, obs group)')
    ap.add_argument('--tenx-compare-csv', action='append', default=[],
                    help='repeatable TEST=path to '
                         'V1_celltyping_compare_Wilcoxon_<TEST>_vs_H2O_veh CSV '
                         '(cols celltype,gene,logfc,method). One per contrast.')
    ap.add_argument('--compare-method', default='marker-typed',
                    help="which typing method to read from the 10X compare CSV "
                         "(default 'marker-typed' so only QC differs)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = list(dict.fromkeys(args.genes))
    norm = Normalize(vmin=args.vmin, vmax=args.vmax, clip=True)
    cmap = plt.get_cmap(args.cmap)
    base_cmap = plt.get_cmap(args.baseline_cmap)
    rng = f'{fmt(args.vmin)}to{fmt(args.vmax)}'

    def baseline_norm_from(df):
        vals = df['mean_expr'].dropna().values
        bvmax = float(np.nanpercentile(vals, 99)) if vals.size else 1.0
        if not np.isfinite(bvmax) or bvmax <= 0:
            bvmax = 1.0
        return Normalize(vmin=0.0, vmax=bvmax, clip=True)

    if args.mode == 'combined':
        req = ['slide_a_dir', 'slide_b_dir', 'slide_a_ann', 'slide_b_ann',
               'wholecell_de_csv', 'loader_script_dir', 'nucleus_h5ad',
               'nucleus_de_csv']
        for r in req:
            if getattr(args, r) is None:
                ap.error(f'--{r.replace("_", "-")} required for --mode combined')
        tables, panel_counts, present_all = [], {}, set()
        loaders = [
            ('WholeCell', load_wholecell,
             lambda: lfc_lookup_wholecell(args.wholecell_de_csv)),
            ('NucleusOnly', load_nucleus,
             lambda: lfc_lookup_nucleus(args.nucleus_de_csv)),
        ]
        for seg, load_fn, lut_fn in loaders:
            logging.info('=== combined: loading %s ===', seg)
            adata = load_fn(args)
            lfc_lut = lut_fn()
            frac_pos, mean_expr, n_cells, present = frac_positive(adata, genes)
            present_all |= present
            for ct in PANELS:
                panel_counts[(seg, ct)] = sum(n_cells.get((ct, g), 0)
                                              for g in GROUPS)
            tables.append(build_table(genes, frac_pos, mean_expr, n_cells,
                                      present, lfc_lut, args.min_cells,
                                      segment=seg))
            del adata
        df = pd.concat(tables, ignore_index=True)
        absent = [g for g in genes if g not in present_all]
        if absent:
            logging.warning('Genes not in either panel: %s', absent)
        logging.info('panel cell counts: %s', panel_counts)
        base = (f'CelltypeIteration_Dotplot_NeuronAstro_Combined_Wilcoxon_'
                f'{args.cmap}_{rng}_{len(genes)}gene-{args.date}')
        df.to_csv(args.outdir / f'{base}.csv', index=False)
        logging.info('Wrote %s.csv', base)
        draw_combined(df, genes, panel_counts, norm, cmap, args.vmin, args.vmax,
                      present_all, args.outdir / f'{base}.png',
                      args.outdir / f'{base}.pdf',
                      baseline_norm_from(df), base_cmap)
        logging.info('Wrote %s.png | .pdf', base)
        logging.info('Done (combined). %d genes (%d absent).',
                     len(genes), len(absent))
        return

    if args.mode == 'pipeline_compare':
        req = ['slide_a_dir', 'slide_b_dir', 'slide_a_ann', 'slide_b_ann',
               'wholecell_de_csv', 'loader_script_dir', 'tenx_h5ad']
        for r in req:
            if getattr(args, r) is None:
                ap.error(f'--{r.replace("_", "-")} required for '
                         '--mode pipeline_compare')
        if not args.tenx_compare_csv:
            ap.error('--tenx-compare-csv required (TEST=path) for '
                     '--mode pipeline_compare')
        tenx_csvs = {}
        for spec in args.tenx_compare_csv:
            if '=' not in spec:
                ap.error(f'--tenx-compare-csv must be TEST=path, got: {spec}')
            t, p = spec.split('=', 1)
            tenx_csvs[t] = Path(p)
        # restrict to contrasts that have a 10X compare table
        contrasts = [c for c in PIPELINE_CONTRASTS if c[0] in tenx_csvs]
        pgroups = sorted(set([t for t, _r, _l in contrasts] + ['H2O_veh']))
        logging.info('pipeline_compare contrasts: %s',
                     [c[2] for c in contrasts])

        tables, panel_counts, present_all = [], {}, set()
        loaders = [
            ('V1_DE', load_wholecell,
             lambda: lfc_lookup_wholecell(args.wholecell_de_csv)),
            ('TenX', load_10x,
             lambda: lfc_lookup_compare(tenx_csvs, args.compare_method)),
        ]
        for seg, load_fn, lut_fn in loaders:
            logging.info('=== pipeline_compare: loading %s ===', seg)
            adata = load_fn(args)
            lfc_lut = lut_fn()
            frac_pos, mean_expr, n_cells, present = frac_positive(adata, genes)
            present_all |= present
            for ct in PANELS:
                panel_counts[(seg, ct)] = sum(n_cells.get((ct, g), 0)
                                              for g in pgroups)
            tables.append(build_table(genes, frac_pos, mean_expr, n_cells,
                                      present, lfc_lut, args.min_cells,
                                      segment=seg, contrasts=contrasts))
            del adata
        df = pd.concat(tables, ignore_index=True)
        absent = [g for g in genes if g not in present_all]
        if absent:
            logging.warning('Genes not in either pipeline: %s', absent)
        logging.info('panel cell counts: %s', panel_counts)
        caption = ('dot size = % cells positive · FC cols colour = Wilcoxon '
                   f'log2FC clipped [{args.vmin:g}, {args.vmax:g}] · baseline '
                   'col colour = mean expr · only difference between halves = '
                   'QC parameters (same marker typing) · † = gene absent '
                   '(dropped by that pipeline’s gene filter)')
        base = (f'CelltypeIteration_Dotplot_NeuronAstro_PipelineCompare_'
                f'V1DEvs10X_Wilcoxon_{args.cmap}_{rng}_{len(genes)}gene-'
                f'{args.date}')
        df.to_csv(args.outdir / f'{base}.csv', index=False)
        logging.info('Wrote %s.csv', base)
        draw_combined(df, genes, panel_counts, norm, cmap, args.vmin, args.vmax,
                      present_all, args.outdir / f'{base}.png',
                      args.outdir / f'{base}.pdf',
                      baseline_norm_from(df), base_cmap,
                      panels=PIPELINE_PANELS, contrasts=contrasts,
                      half_specs=HALF_SPECS_PIPE, caption=caption)
        logging.info('Wrote %s.png | .pdf', base)
        logging.info('Done (pipeline_compare). %d genes (%d absent).',
                     len(genes), len(absent))
        return

    if args.mode == 'wholecell':
        for req in ['slide_a_dir', 'slide_b_dir', 'slide_a_ann', 'slide_b_ann',
                    'wholecell_de_csv', 'loader_script_dir']:
            if getattr(args, req) is None:
                ap.error(f'--{req.replace("_", "-")} required for --mode wholecell')
        logging.info('=== whole-cell: load + marker celltyping (pool A+B) ===')
        adata = load_wholecell(args)
        lfc_lut = lfc_lookup_wholecell(args.wholecell_de_csv)
        title = 'Neuron vs Astrocyte — WHOLE-CELL iteration (dot plot)'
        seg = 'WholeCell'
    else:
        for req in ['nucleus_h5ad', 'nucleus_de_csv']:
            if getattr(args, req) is None:
                ap.error(f'--{req.replace("_", "-")} required for --mode nucleus')
        logging.info('=== nucleus: load h5ad ===')
        adata = load_nucleus(args)
        lfc_lut = lfc_lookup_nucleus(args.nucleus_de_csv)
        title = 'Neuron vs Astrocyte — NUCLEUS-ONLY iteration (dot plot)'
        seg = 'NucleusOnly'

    frac_pos, mean_expr, n_cells, present = frac_positive(adata, genes)
    absent = [g for g in genes if g not in present]
    if absent:
        logging.warning('Genes not in panel: %s', absent)
    panel_counts = {ct: sum(n_cells.get((ct, g), 0) for g in GROUPS) for ct in PANELS}
    logging.info('panel cell counts (contrast groups): %s', panel_counts)

    df = build_table(genes, frac_pos, mean_expr, n_cells, present, lfc_lut,
                     args.min_cells, segment=seg)

    base = (f'CelltypeIteration_Dotplot_NeuronAstro_{seg}_Wilcoxon_'
            f'{args.cmap}_{rng}_{len(genes)}gene-{args.date}')
    csv_path = args.outdir / f'{base}.csv'
    df.to_csv(csv_path, index=False)
    logging.info('Wrote %s', csv_path.name)

    out_png = args.outdir / f'{base}.png'
    out_pdf = args.outdir / f'{base}.pdf'
    draw(df, genes, panel_counts, norm, cmap, args.vmin, args.vmax, present,
         title, out_png, out_pdf, baseline_norm_from(df), base_cmap)
    logging.info('Wrote %s | %s', out_png.name, out_pdf.name)
    logging.info('Done. %d genes (%d absent), colourbar [%g, %g].',
                 len(genes), len(absent), args.vmin, args.vmax)


if __name__ == '__main__':
    main()
