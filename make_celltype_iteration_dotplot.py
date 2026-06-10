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
    ('Neuronal markers',   ['Rbfox3', 'Snap25', 'Syn1', 'Stmn2', 'Map2', 'Tubb3']),
    ('Astrocyte markers',  ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1']),
    ('Activity / receptors', ['Grin2d', 'Oprm1', 'Bdnf', 'Arc', 'Egr1', 'Fos']),
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
    and n cells per (cell_type, group). Uses layers['counts'] if present."""
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
    frac_pos = {}
    for gene in present:
        c = dense_col(sub[:, gene].layers['counts'] if use_counts else sub[:, gene].X)
        pos = (c > 0).astype(float)
        for k, idx in keyidx.items():
            if idx.size:
                frac_pos[(gene, *k)] = float(pos[idx].mean())
    return frac_pos, n_cells, set(present)


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
def build_table(genes, frac_pos, n_cells, present, lfc_lut, min_cells):
    rows = []
    for gene in genes:
        in_panel = gene in present
        for ct in PANELS:
            for test, ref, lab in CONTRASTS:
                n_t = n_cells.get((ct, test), 0)
                fv = frac_pos.get((gene, ct, test), np.nan)
                lfc = lfc_lut.get((gene, ct, test), np.nan)
                ok = in_panel and n_t >= min_cells
                rows.append({
                    'gene': gene, 'in_panel': in_panel, 'cell_type': ct,
                    'contrast': lab, 'test': test, 'reference': ref,
                    'n_test_cells': int(n_t),
                    'frac_pos_test': (float(fv) if ok and not pd.isna(fv) else np.nan),
                    'pct_pos_test': (float(fv) * 100.0 if ok and not pd.isna(fv) else np.nan),
                    'log2FC': (float(lfc) if ok and not pd.isna(lfc) else np.nan),
                })
    return pd.DataFrame(rows)


# --------------------------- drawing ----------------------------------------
def draw(df, genes, panel_counts, norm, cmap, vmin, vmax, present, title,
         out_png, out_pdf, gene_fs=8.5, col_fs=8.0, panel_fs=11.0,
         size_scale=1.0):
    n = len(genes)
    y_of = {g: i for i, g in enumerate(genes)}
    ncol = len(CONTRASTS)
    col_labels = [lab for _t, _r, lab in CONTRASTS]

    def sz(frac):
        return size_of(frac) * size_scale

    fig_w = 2.6 + len(PANELS) * (0.52 * ncol + 0.6)
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
        for gene, yi in y_of.items():
            for j, (_t, _r, lab) in enumerate(CONTRASTS):
                fv, lfc = lut.get((gene, lab), (np.nan, np.nan))
                if pd.isna(fv) or pd.isna(lfc) or fv <= 0:
                    continue
                ax.scatter(j, yi, s=sz(fv), c=[lfc], cmap=cmap, norm=norm,
                           edgecolor='black', linewidth=0.4, zorder=3)
        ax.set_xlim(-0.6, ncol - 0.4)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks(range(ncol))
        ax.set_xticklabels(col_labels, rotation=45, ha='right', va='top', fontsize=col_fs)
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

    # category labels down the far left
    block_start = 0
    for lab, gl in GENE_GROUPS:
        present_in_group = [g for g in gl if g in genes]
        if not present_in_group:
            continue
        mid = block_start + (len(present_in_group) - 1) / 2.0
        axes[0].text(-1.75, mid, lab, rotation=90, ha='center', va='center',
                     fontsize=7.5, fontweight='bold', color='0.35',
                     transform=axes[0].transData)
        block_start += len(present_in_group)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.024, pad=0.04, extend='both')
    cbar.set_label('Wilcoxon log2 fold change (test / H2O_veh)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    handles = [ax.scatter([], [], s=sz(f), c='0.55', edgecolor='black',
                          linewidth=0.4) for f in LEGEND_FRACS]
    labels = [f'{int(f * 100)}%' for f in LEGEND_FRACS]
    fig.legend(handles, labels, title='% cells\npositive',
               loc='center left', bbox_to_anchor=(0.995, 0.30), frameon=False,
               labelspacing=1.6, handletextpad=1.0, borderpad=0.8,
               fontsize=8, title_fontsize=8)

    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.998)
    fig.text(0.5, 0.972,
             'dot size = % cells positive (test group) · colour = Wilcoxon log2FC, '
             f'clipped [{vmin:g}, {vmax:g}]  ·  † = gene absent from panel',
             ha='center', va='top', fontsize=7.5, color='0.35')

    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


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


# --------------------------- main -------------------------------------------
def fmt(v):
    return ('%g' % v).replace('-', 'neg').replace('.', 'p')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', required=True, choices=['wholecell', 'nucleus'])
    ap.add_argument('--genes', nargs='+', default=DEFAULT_GENES)
    ap.add_argument('--vmin', type=float, default=-1.0)
    ap.add_argument('--vmax', type=float, default=1.0)
    ap.add_argument('--cmap', default='RdBu_r')
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
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = list(dict.fromkeys(args.genes))

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

    frac_pos, n_cells, present = frac_positive(adata, genes)
    absent = [g for g in genes if g not in present]
    if absent:
        logging.warning('Genes not in panel: %s', absent)
    panel_counts = {ct: sum(n_cells.get((ct, g), 0) for g in GROUPS) for ct in PANELS}
    logging.info('panel cell counts (contrast groups): %s', panel_counts)

    df = build_table(genes, frac_pos, n_cells, present, lfc_lut, args.min_cells)

    rng = f'{fmt(args.vmin)}to{fmt(args.vmax)}'
    base = (f'CelltypeIteration_Dotplot_NeuronAstro_{seg}_Wilcoxon_'
            f'{args.cmap}_{rng}_{len(genes)}gene-{args.date}')
    csv_path = args.outdir / f'{base}.csv'
    df.to_csv(csv_path, index=False)
    logging.info('Wrote %s', csv_path.name)

    norm = Normalize(vmin=args.vmin, vmax=args.vmax, clip=True)
    cmap = plt.get_cmap(args.cmap)
    out_png = args.outdir / f'{base}.png'
    out_pdf = args.outdir / f'{base}.pdf'
    draw(df, genes, panel_counts, norm, cmap, args.vmin, args.vmax, present,
         title, out_png, out_pdf)
    logging.info('Wrote %s | %s', out_png.name, out_pdf.name)
    logging.info('Done. %d genes (%d absent), colourbar [%g, %g].',
                 len(genes), len(absent), args.vmin, args.vmax)


if __name__ == '__main__':
    main()
