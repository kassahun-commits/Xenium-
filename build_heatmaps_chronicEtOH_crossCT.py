#!/usr/bin/env python3
"""
Cross-cell-type heat maps for the CHRONIC ethanol contrast, Xenium May 2026
hippocampus dataset.

Single contrast (one column per cell-type block):
    Chronic EtOH = ChronicEtOH vs H2O_veh

Each page anchors on one cell type's significant-gene list, then shows those
SAME genes' log2FC across EVERY cell type (one column per cell type, placed
right next to each other) so each gene can be read straight across cell types.

Gene (row) selection (in the anchor cell type):
    significant in ChronicEtOH vs H2O_veh (|log2FC| >= LFC_THRESH AND padj < PADJ_THRESH),
    then top TOP_N up + TOP_N down by log2FC.

Variants:
    3CT  : cell types Excitatory, Inhibitory, Astrocyte ; anchors Excitatory, Astrocyte
           (Inhibitory anchor skipped: ~0 significant genes)
    NvsA : cell types Neuron, Astrocyte                 ; anchors Neuron, Astrocyte
Shown for V1 and V2.
Output : one multi-page summary PDF + per-page source-data CSVs.
"""
from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

LFC_THRESH = 0.5
PADJ_THRESH = 1e-3
CLIP = 3.0
TOP_N = 50

CONTRAST = ('ChronicEtOH', 'H2O_veh', 'Chronic EtOH')
SEL_TEST, SEL_REF = CONTRAST[0], CONTRAST[1]

# variant -> (cell-type blocks shown, anchor cell types)
VARIANTS = {
    '3CT':  (['Excitatory', 'Inhibitory', 'Astrocyte'], ['Excitatory', 'Astrocyte']),
    'NvsA': (['Neuron', 'Astrocyte'], ['Neuron', 'Astrocyte']),
}


def cl(df, subtype, test, ref):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['logfc']


def cp(df, subtype, test, ref):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['padj']


def anchor_genes(df, anchor):
    lfc = cl(df, anchor, SEL_TEST, SEL_REF)
    padj = cp(df, anchor, SEL_TEST, SEL_REF)
    if len(lfc) == 0:
        return []
    sig = (lfc.abs() >= LFC_THRESH) & (padj < PADJ_THRESH)
    sel = lfc[sig]
    up = sel[sel > 0].sort_values(ascending=False).head(TOP_N)
    dn = sel[sel < 0].sort_values(ascending=True).head(TOP_N)
    return list(up.index) + list(dn.index[::-1])


def build_matrix(df, genes, celltypes):
    """rows=genes, one column per cell type (Chronic EtOH log2FC)."""
    cols = {ct: cl(df, ct, SEL_TEST, SEL_REF).reindex(genes) for ct in celltypes}
    return pd.DataFrame(cols, index=genes)


def draw_page(pdf, df, version, variant, anchor, csv_path):
    celltypes, _ = VARIANTS[variant]
    genes = anchor_genes(df, anchor)
    nct = len(celltypes)
    fig_w = 2.2 + 0.75 * nct
    fig_h = max(6.0, 0.135 * max(len(genes), 1) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    title = (f'{version} — Chronic EtOH vs H2O_veh — anchored on {anchor} genes\n'
             f'rows = top {TOP_N} up + {TOP_N} down by log2FC in {anchor}; '
             f'one column per cell type; log2FC clipped ±{CLIP:g}')

    if not genes:
        ax.text(0.5, 0.5, f'(no genes for {anchor})', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(title, fontsize=9, y=0.99)
        pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)
        return

    mat = build_matrix(df, genes, celltypes)
    vals = mat.values.astype(float)
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad('lightgray')
    im = ax.imshow(np.ma.masked_invalid(vals), aspect='auto', cmap=cmap,
                   vmin=-CLIP, vmax=CLIP)

    ax.set_xticks(range(nct))
    ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=8, fontweight='bold')
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=max(3, min(6, 520 / max(len(genes), 1))))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # separators between cell-type columns
    for b in range(1, nct):
        ax.axvline(b - 0.5, color='white', lw=2.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label('log2 fold change', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=9, y=1.0)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    out = mat.copy()
    out.columns = [f'{ct}|Chronic EtOH' for ct in mat.columns]
    out.insert(0, 'anchor_gene', out.index)
    out.to_csv(csv_path, index=False)


def cover(title, body, path):
    d = fitz.open(); p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 60, 738, 160), title, fontsize=17,
                     fontname='Helvetica-Bold', align=1)
    p.insert_textbox(fitz.Rect(54, 180, 738, 560), body, fontsize=11,
                     fontname='Helvetica')
    d.save(str(path)); d.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--nvsa-v1', required=True, type=Path)
    ap.add_argument('--nvsa-v2', required=True, type=Path)
    ap.add_argument('--ct3-v1', required=True, type=Path)
    ap.add_argument('--ct3-v2', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = {
        ('NvsA', 'V1'): pd.read_csv(args.nvsa_v1),
        ('NvsA', 'V2'): pd.read_csv(args.nvsa_v2),
        ('3CT', 'V1'): pd.read_csv(args.ct3_v1),
        ('3CT', 'V2'): pd.read_csv(args.ct3_v2),
    }

    tmp = args.outdir / '_tmp_hmx_chronic'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    def jobs():
        for version in ('V1', 'V2'):
            for variant in ('3CT', 'NvsA'):
                _, anchors = VARIANTS[variant]
                for anchor in anchors:
                    yield version, variant, anchor

    n = 0
    with PdfPages(pages_pdf) as pdf:
        for version, variant, anchor in jobs():
            csv_path = (args.outdir /
                        f'Heatmap_chronicEtOH_crossCT_{version}_{variant}_anchor{anchor}_{today}.csv')
            draw_page(pdf, data[(variant, version)], version, variant, anchor, csv_path)
            logging.info('Page: %s %s anchor=%s', version, variant, anchor)
            n += 1

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — Cross-cell-type heat maps\n(Chronic EtOH vs H2O_veh)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "Each page anchors on one cell type's significant-gene list, then shows those\n"
          "SAME genes' log2FC across every cell type (one column per cell type, side by\n"
          "side) so each gene can be read straight across cell types.\n\n"
          "Single contrast: Chronic EtOH (ChronicEtOH) vs H2O_veh.\n"
          f"Cell value = log2 fold change; diverging blue-white-red, clipped to ±{CLIP:g}.\n"
          f"Row selection: |log2FC| >= {LFC_THRESH:g} AND adj p < {PADJ_THRESH:g} in the anchor\n"
          f"cell type; top {TOP_N} up + {TOP_N} down by log2FC.\n\n"
          "Variants: 3CT (columns Excitatory, Inhibitory, Astrocyte; anchors Excitatory &\n"
          "Astrocyte; Inhibitory anchor omitted — ~0 significant genes) and NvsA (columns\n"
          "Neuron, Astrocyte; anchors Neuron & Astrocyte). Shown for V1 and V2.\n"
          "Gray cells = gene not present / not tested in that cell type.\n\n"
          "Statistics caveat: cell-level Wilcoxon DE — adj p-values inflated by\n"
          "pseudoreplication; treat direction and log2FC magnitude as primary.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=n - 1)
    out_pdf = args.outdir / f'Xenium_Heatmaps_chronicEtOH_crossCellType_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + n)

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
