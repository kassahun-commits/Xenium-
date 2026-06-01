#!/usr/bin/env python3
"""
Heat maps of log2FC for the CHRONIC ethanol contrast, Xenium May 2026
hippocampus dataset.

Single contrast (one column per cell-type panel):
    Chronic EtOH = ChronicEtOH vs H2O_veh

Row (gene) selection (done independently within each cell type):
    genes significant in ChronicEtOH vs H2O_veh
    (|log2FC| >= LFC_THRESH AND padj < PADJ_THRESH),
    then top TOP_N up- and top TOP_N down-regulated by that log2FC.

Variants : all cell types (Excit/Inhib/Astro) and Neuron+Astrocyte, each V1 & V2.
Output   : one multi-page summary PDF + per-page source-data CSVs.
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
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

LFC_THRESH = 0.5
PADJ_THRESH = 1e-3
CLIP = 3.0
TOP_N = 50

# (test, reference, column label)
CONTRAST = ('ChronicEtOH', 'H2O_veh', 'Chronic EtOH')
SEL_TEST, SEL_REF = CONTRAST[0], CONTRAST[1]

VARIANTS = {
    '3CT':  ['Excitatory', 'Inhibitory', 'Astrocyte'],
    'NvsA': ['Neuron', 'Astrocyte'],
}


def contrast_series(df, subtype, test, ref, col):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')[col]


def select_genes(df, subtype):
    """Return ordered list of genes (top up then top down) and the 1-col matrix."""
    lfc = contrast_series(df, subtype, SEL_TEST, SEL_REF, 'logfc')
    padj = contrast_series(df, subtype, SEL_TEST, SEL_REF, 'padj')
    if len(lfc) == 0:
        return [], None
    sig = (lfc.abs() >= LFC_THRESH) & (padj < PADJ_THRESH)
    sel = lfc[sig]
    up = sel[sel > 0].sort_values(ascending=False).head(TOP_N)
    dn = sel[sel < 0].sort_values(ascending=True).head(TOP_N)
    ordered = list(up.index) + list(dn.index[::-1])  # strongest up top, strongest down bottom
    mat = pd.DataFrame({CONTRAST[2]: lfc.reindex(ordered)}, index=ordered)
    return ordered, mat


def draw_page(pdf, df, version, variant, csv_path):
    subtypes = VARIANTS[variant]
    panels = [(st, *select_genes(df, st)) for st in subtypes]
    maxg = max((len(g) for _, g, _ in panels), default=1)
    n = len(subtypes)
    fig_w = 1.9 * n + 1.6
    fig_h = max(6.0, 0.135 * maxg + 2.2)
    fig, axes = plt.subplots(1, n, figsize=(fig_w, fig_h))
    if n == 1:
        axes = [axes]
    im = None
    rows_out = []
    for ax, (st, genes, mat) in zip(axes, panels):
        ax.set_title(f'{st}\n(n={len(genes)})', fontsize=9, fontweight='bold')
        if not genes:
            ax.text(0.5, 0.5, '(no genes)', ha='center', va='center',
                    transform=ax.transAxes, fontsize=8, color='gray')
            ax.set_xticks([]); ax.set_yticks([])
            continue
        vals = mat.values.astype(float)
        im = ax.imshow(vals, aspect='auto', cmap='RdBu_r', vmin=-CLIP, vmax=CLIP)
        ax.set_xticks(range(len(mat.columns)))
        ax.set_xticklabels(mat.columns, rotation=45, ha='right', fontsize=7)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=max(3, min(6, 520 / max(maxg, 1))))
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        m2 = mat.copy()
        m2.insert(0, 'gene', m2.index)
        m2.insert(0, 'cell_type', st)
        rows_out.append(m2)
    var_name = ('All cell types (Excit/Inhib/Astro)' if variant == '3CT'
                else 'Neuron vs Astrocyte')
    fig.suptitle(f'{version} — {var_name}\n'
                 f'Chronic EtOH vs H2O_veh log2FC; rows = top {TOP_N} up + top {TOP_N} down '
                 f'by log2FC (per cell type); color clipped ±{CLIP:g}',
                 fontsize=10, y=0.995)
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.03)
        cbar.set_label('log2 fold change', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    plt.tight_layout(rect=[0, 0, 0.9, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)
    if rows_out:
        pd.concat(rows_out).to_csv(csv_path, index=False)


def cover(text_title, body, path):
    d = fitz.open()
    p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 70, 738, 150), text_title,
                     fontsize=18, fontname='Helvetica-Bold', align=1)
    p.insert_textbox(fitz.Rect(54, 170, 738, 560), body,
                     fontsize=11, fontname='Helvetica')
    d.save(str(path)); d.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--nvsa-v1', required=True, type=Path)
    ap.add_argument('--nvsa-v2', required=True, type=Path)
    ap.add_argument('--ct3-v1', required=True, type=Path)
    ap.add_argument('--ct3-v2', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = {
        ('NvsA', 'V1'): pd.read_csv(args.nvsa_v1),
        ('NvsA', 'V2'): pd.read_csv(args.nvsa_v2),
        ('3CT', 'V1'): pd.read_csv(args.ct3_v1),
        ('3CT', 'V2'): pd.read_csv(args.ct3_v2),
    }

    tmp = args.outdir / '_tmp_hm_chronic'
    tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    order = [('V1', '3CT'), ('V1', 'NvsA'), ('V2', '3CT'), ('V2', 'NvsA')]
    with PdfPages(pages_pdf) as pdf:
        for version, variant in order:
            csv_path = (args.outdir /
                        f'Heatmap_chronicEtOH_{version}_{variant}_{today}.csv')
            draw_page(pdf, data[(variant, version)], version, variant, csv_path)
            logging.info('Page done: %s %s -> %s', version, variant, csv_path.name)

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — Heat maps (Chronic EtOH vs H2O_veh)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          f"Single contrast: Chronic EtOH (ChronicEtOH) vs H2O_veh.\n"
          f"Cell value = log2 fold change; diverging blue-white-red, clipped to ±{CLIP:g}.\n\n"
          f"Gene (row) selection (independently within each cell type):\n"
          f"  significant in Chronic EtOH vs H2O_veh: |log2FC| >= {LFC_THRESH:g} AND adj p < {PADJ_THRESH:g}.\n"
          f"  Rows = top {TOP_N} up + top {TOP_N} down ranked by log2FC.\n\n"
          "Two variants per version: all cell types (Excitatory/Inhibitory/Astrocyte),\n"
          "and Neuron vs Astrocyte. Shown for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).\n\n"
          "Statistics caveat: cell-level Wilcoxon DE — adj p-values inflated by\n"
          "pseudoreplication; treat direction and log2FC magnitude as primary.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=3)
    out_pdf = args.outdir / f'Xenium_Heatmaps_chronicEtOH_Summary_{today}.pdf'
    out.save(str(out_pdf))
    out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + 4)

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
