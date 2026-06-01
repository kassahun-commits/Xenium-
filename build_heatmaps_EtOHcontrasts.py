#!/usr/bin/env python3
"""
Heat maps of log2FC across three H2O-referenced contrasts, for the Xenium
May 2026 hippocampus dataset.

Columns of every heat map (all vs H2O_veh):
    EtOH         = EtOH_veh    vs H2O_veh
    MCT1i        = H2O_MCT1i   vs H2O_veh
    EtOH+MCT1i   = EtOH_MCT1i  vs H2O_veh

Row (gene) selection:
    Set 1  "EtOH-significant"  : genes significant in EtOH vs H2O
                                 (|log2FC| >= LFC_THRESH AND padj < PADJ_THRESH)
    Set 2  "EtOH-specific"     : Set-1 genes that do NOT also pass that same
                                 significance bar in MCT1i vs H2O
Within each cell type we keep the top TOP_N up- and top TOP_N down-regulated
genes (ranked by EtOH log2FC).

Variants : Neuron+Astrocyte (NvsA) and Excit+Inhib+Astro (3CT), each for V1 & V2.
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
import numpy as np
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

LFC_THRESH = 0.5
PADJ_THRESH = 1e-3
CLIP = 3.0
TOP_N = 50

# (test, reference, column label)
CONTRASTS = [
    ('EtOH_veh',   'H2O_veh', 'EtOH'),
    ('H2O_MCT1i',  'H2O_veh', 'MCT1i'),
    ('EtOH_MCT1i', 'H2O_veh', 'EtOH+MCT1i'),
]
SEL_TEST, SEL_REF = 'EtOH_veh', 'H2O_veh'      # selection contrast
MCT_TEST, MCT_REF = 'H2O_MCT1i', 'H2O_veh'     # "did not change" contrast (Set 2)

VARIANTS = {
    'NvsA': ['Neuron', 'Astrocyte'],
    '3CT':  ['Excitatory', 'Inhibitory', 'Astrocyte'],
}


def contrast_lfc(df, subtype, test, ref):
    """Return Series gene -> log2FC for one subtype/contrast."""
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['logfc']


def contrast_padj(df, subtype, test, ref):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['padj']


def select_genes(df, subtype, which_set):
    """Return ordered list of genes (top up then top down) and the value matrix."""
    e_lfc = contrast_lfc(df, subtype, SEL_TEST, SEL_REF)
    e_padj = contrast_padj(df, subtype, SEL_TEST, SEL_REF)
    if len(e_lfc) == 0:
        return [], None
    sig = (e_lfc.abs() >= LFC_THRESH) & (e_padj < PADJ_THRESH)
    genes = e_lfc.index[sig]
    if which_set == 'set2':
        m_lfc = contrast_lfc(df, subtype, MCT_TEST, MCT_REF)
        m_padj = contrast_padj(df, subtype, MCT_TEST, MCT_REF)
        m_sig = (m_lfc.abs() >= LFC_THRESH) & (m_padj < PADJ_THRESH)
        m_sig = m_sig.reindex(genes).fillna(False)
        genes = genes[~m_sig.values]
    sel = e_lfc.loc[genes]
    up = sel[sel > 0].sort_values(ascending=False).head(TOP_N)
    dn = sel[sel < 0].sort_values(ascending=True).head(TOP_N)
    ordered = list(up.index) + list(dn.index[::-1])  # strongest up top, strongest down bottom
    # build value matrix (genes x contrasts) of log2FC
    cols = {}
    for test, ref, lab in CONTRASTS:
        cols[lab] = contrast_lfc(df, subtype, test, ref).reindex(ordered)
    mat = pd.DataFrame(cols, index=ordered)
    return ordered, mat


def draw_page(pdf, dfs_for_version, version, variant, which_set, csv_path):
    subtypes = VARIANTS[variant]
    df = dfs_for_version
    panels = [(st, *select_genes(df, st, which_set)) for st in subtypes]
    maxg = max((len(g) for _, g, _ in panels), default=1)
    n = len(subtypes)
    fig_w = 2.6 * n + 1.4
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
    set_name = ('EtOH-significant genes' if which_set == 'set1'
                else 'EtOH-specific (unchanged by MCT1i)')
    var_name = 'Neuron vs Astrocyte' if variant == 'NvsA' else 'All cell types (Excit/Inhib/Astro)'
    fig.suptitle(f'{version} — {var_name}\nSet: {set_name}   |   '
                 f'log2FC vs H2O_veh; rows = top {TOP_N} up + top {TOP_N} down by EtOH log2FC; '
                 f'color clipped ±{CLIP:g}',
                 fontsize=10, y=0.995)
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label('log2 fold change', fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    plt.tight_layout(rect=[0, 0, 0.93, 0.96])
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


def divider(title, path):
    d = fitz.open()
    p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 250, 738, 360), title,
                     fontsize=20, fontname='Helvetica-Bold', align=1)
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

    tmp = args.outdir / '_tmp_hm'
    tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    # page order: per set -> V1 3CT, V1 NvsA, V2 3CT, V2 NvsA
    order = [('V1', '3CT'), ('V1', 'NvsA'), ('V2', '3CT'), ('V2', 'NvsA')]
    with PdfPages(pages_pdf) as pdf:
        for which_set in ('set1', 'set2'):
            for version, variant in order:
                csv_path = (args.outdir /
                            f'Heatmap_{which_set}_{version}_{variant}_{today}.csv')
                draw_page(pdf, data[(variant, version)], version, variant,
                          which_set, csv_path)
                logging.info('Page done: %s %s %s -> %s',
                             which_set, version, variant, csv_path.name)

    # assemble final PDF with cover + dividers
    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — Heat maps (EtOH / MCT1i / EtOH+MCT1i vs H2O)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          f"Columns (all vs H2O_veh): EtOH, MCT1i, EtOH+MCT1i.\n"
          f"Cell value = log2 fold change; diverging blue-white-red, clipped to ±{CLIP:g}.\n\n"
          f"Gene (row) selection thresholds: |log2FC| >= {LFC_THRESH:g} AND adj p < {PADJ_THRESH:g}.\n"
          f"Rows = top {TOP_N} up + top {TOP_N} down ranked by EtOH log2FC (per cell type).\n\n"
          "Set 1 — EtOH-significant: genes significant in EtOH vs H2O.\n"
          "Set 2 — EtOH-specific: Set-1 genes that are NOT significant in MCT1i vs H2O\n"
          "        (i.e. changed by ethanol but not by drug alone).\n\n"
          "Variants: Neuron vs Astrocyte, and all cell types (Excitatory/Inhibitory/\n"
          "Astrocyte). Shown for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).\n\n"
          "Statistics caveat: cell-level Wilcoxon DE — adj p-values inflated by\n"
          "pseudoreplication; treat direction and log2FC magnitude as primary.",
          cover_path)
    div1 = tmp / 'div1.pdf'; divider('SET 1 — EtOH-significant genes', div1)
    div2 = tmp / 'div2.pdf'; divider('SET 2 — EtOH-specific genes\n(unchanged by MCT1i)', div2)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(fitz.open(str(div1)))
    out.insert_pdf(src, from_page=0, to_page=3)   # set1: 4 pages
    out.insert_pdf(fitz.open(str(div2)))
    out.insert_pdf(src, from_page=4, to_page=7)   # set2: 4 pages
    out_pdf = args.outdir / f'Xenium_Heatmaps_EtOHcontrasts_Summary_{today}.pdf'
    out.save(str(out_pdf))
    out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + 1 + 4 + 1 + 4)

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
