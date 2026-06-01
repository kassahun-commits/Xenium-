#!/usr/bin/env python3
"""
Cross-cell-type heat maps: pick a gene set from one "anchor" cell type, then
show how those SAME genes change across EVERY cell type, over four contrasts.

Columns (per cell-type block):
    EtOH            = EtOH_veh    vs H2O_veh
    MCT1i           = H2O_MCT1i   vs H2O_veh
    EtOH+MCT1i      = EtOH_MCT1i  vs H2O_veh
    EtOH+MCT1i/EtOH = EtOH_MCT1i  vs EtOH_veh   (direct MCT1i rescue contrast)

Each page = one heat map: rows = anchor cell type's top TOP_N up + TOP_N down
genes (ranked by EtOH log2FC), columns grouped into one block per cell type so
each gene can be read straight across all cell types.

Gene (row) selection thresholds: |log2FC| >= LFC_THRESH AND padj < PADJ_THRESH
in the anchor cell type's EtOH-vs-H2O contrast.
    Set 1 "EtOH-significant" : anchor genes significant in EtOH vs H2O.
    Set 2 "EtOH-specific"    : Set-1 genes NOT significant in MCT1i vs H2O
                               (in the anchor cell type).

Variants:
    NvsA : cell types Neuron, Astrocyte          ; anchors Neuron, Astrocyte
    3CT  : cell types Excitatory, Inhibitory,    ; anchors Excitatory, Astrocyte
           Astrocyte                               (Inhibitory anchor skipped: ~0 genes)
Shown for V1 and V2.
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

CONTRASTS = [
    ('EtOH_veh',   'H2O_veh',  'EtOH'),
    ('H2O_MCT1i',  'H2O_veh',  'MCT1i'),
    ('EtOH_MCT1i', 'H2O_veh',  'EtOH+MCT1i'),
    ('EtOH_MCT1i', 'EtOH_veh', 'E+M / EtOH'),
]
SEL_TEST, SEL_REF = 'EtOH_veh', 'H2O_veh'
MCT_TEST, MCT_REF = 'H2O_MCT1i', 'H2O_veh'

# variant -> (cell-type blocks shown, anchor cell types)
VARIANTS = {
    'NvsA': (['Neuron', 'Astrocyte'], ['Neuron', 'Astrocyte']),
    '3CT':  (['Excitatory', 'Inhibitory', 'Astrocyte'], ['Excitatory', 'Astrocyte']),
}


def cl(df, subtype, test, ref):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['logfc']


def cp(df, subtype, test, ref):
    s = df[(df['subtype'] == subtype) & (df['test'] == test) & (df['reference'] == ref)]
    return s.set_index('gene')['padj']


def anchor_genes(df, anchor, which_set):
    e_lfc = cl(df, anchor, SEL_TEST, SEL_REF)
    e_padj = cp(df, anchor, SEL_TEST, SEL_REF)
    if len(e_lfc) == 0:
        return []
    sig = (e_lfc.abs() >= LFC_THRESH) & (e_padj < PADJ_THRESH)
    genes = e_lfc.index[sig]
    if which_set == 'set2':
        m_lfc = cl(df, anchor, MCT_TEST, MCT_REF)
        m_padj = cp(df, anchor, MCT_TEST, MCT_REF)
        m_sig = ((m_lfc.abs() >= LFC_THRESH) & (m_padj < PADJ_THRESH)).reindex(genes).fillna(False)
        genes = genes[~m_sig.values]
    sel = e_lfc.loc[genes]
    up = sel[sel > 0].sort_values(ascending=False).head(TOP_N)
    dn = sel[sel < 0].sort_values(ascending=True).head(TOP_N)
    return list(up.index) + list(dn.index[::-1])


def build_matrix(df, genes, celltypes):
    """Return DataFrame rows=genes, MultiIndex cols (celltype, contrast label)."""
    cols = {}
    for ct in celltypes:
        for test, ref, lab in CONTRASTS:
            cols[(ct, lab)] = cl(df, ct, test, ref).reindex(genes)
    mat = pd.DataFrame(cols, index=genes)
    mat.columns = pd.MultiIndex.from_tuples(mat.columns, names=['cell_type', 'contrast'])
    return mat


def draw_page(pdf, df, version, variant, anchor, which_set, csv_path):
    celltypes, _ = VARIANTS[variant]
    genes = anchor_genes(df, anchor, which_set)
    nct = len(celltypes)
    ncon = len(CONTRASTS)
    ncol = nct * ncon
    fig_w = 1.9 + 0.46 * ncol
    fig_h = max(6.0, 0.135 * max(len(genes), 1) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    set_name = ('EtOH-significant' if which_set == 'set1'
                else 'EtOH-specific (unchanged by MCT1i)')
    title = (f'{version} — anchored on {anchor} genes — Set: {set_name}\n'
             f'rows = top {TOP_N} up + {TOP_N} down by EtOH log2FC in {anchor}; '
             f'columns grouped by cell type; log2FC clipped ±{CLIP:g}')

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

    ax.set_xticks(range(ncol))
    ax.set_xticklabels([lab for _ in celltypes for *_, lab in CONTRASTS],
                       rotation=90, fontsize=6)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=max(3, min(6, 520 / max(len(genes), 1))))
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)

    # cell-type block separators + headers
    for b in range(1, nct):
        ax.axvline(b * ncon - 0.5, color='white', lw=2.5)
    for b, ct in enumerate(celltypes):
        ax.text(b * ncon + ncon / 2 - 0.5, 1.012, ct,
                transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                fontsize=8.5, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('log2 fold change', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.suptitle(title, fontsize=9, y=1.0)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    out = mat.copy()
    out.columns = [f'{ct}|{lab}' for ct, lab in mat.columns]
    out.insert(0, 'anchor_gene', out.index)
    out.to_csv(csv_path, index=False)


def cover(title, body, path):
    d = fitz.open(); p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 60, 738, 150), title, fontsize=17,
                     fontname='Helvetica-Bold', align=1)
    p.insert_textbox(fitz.Rect(54, 170, 738, 560), body, fontsize=11,
                     fontname='Helvetica')
    d.save(str(path)); d.close()


def divider(title, path):
    d = fitz.open(); p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 250, 738, 360), title, fontsize=20,
                     fontname='Helvetica-Bold', align=1)
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

    tmp = args.outdir / '_tmp_hmx'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    # collect page jobs in order, per set
    def jobs():
        for version in ('V1', 'V2'):
            for variant in ('3CT', 'NvsA'):
                _, anchors = VARIANTS[variant]
                for anchor in anchors:
                    yield version, variant, anchor

    page_counts = {}
    with PdfPages(pages_pdf) as pdf:
        for which_set in ('set1', 'set2'):
            n = 0
            for version, variant, anchor in jobs():
                csv_path = (args.outdir /
                            f'Heatmap_crossCT_{which_set}_{version}_{variant}_anchor{anchor}_{today}.csv')
                draw_page(pdf, data[(variant, version)], version, variant, anchor,
                          which_set, csv_path)
                logging.info('Page: %s %s %s anchor=%s', which_set, version, variant, anchor)
                n += 1
            page_counts[which_set] = n

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — Cross-cell-type heat maps\n(EtOH / MCT1i / EtOH+MCT1i vs H2O; + EtOH+MCT1i vs EtOH)',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "Each page anchors on one cell type's significant-gene list, then shows those\n"
          "SAME genes across every cell type so each gene can be read across cell types.\n\n"
          "Columns per cell-type block (4): EtOH (vs H2O), MCT1i (vs H2O),\n"
          "EtOH+MCT1i (vs H2O), and E+M/EtOH = EtOH+MCT1i vs EtOH_veh (direct rescue).\n\n"
          f"Cell value = log2 fold change; diverging blue-white-red, clipped to ±{CLIP:g}.\n"
          f"Row selection: |log2FC| >= {LFC_THRESH:g} AND adj p < {PADJ_THRESH:g} in the anchor\n"
          f"cell type's EtOH-vs-H2O contrast; top {TOP_N} up + {TOP_N} down by EtOH log2FC.\n\n"
          "Set 1 — EtOH-significant: anchor genes significant in EtOH vs H2O.\n"
          "Set 2 — EtOH-specific: Set-1 genes NOT significant in MCT1i vs H2O.\n\n"
          "Variants: NvsA (blocks Neuron, Astrocyte; anchors Neuron & Astrocyte) and\n"
          "3CT (blocks Excitatory, Inhibitory, Astrocyte; anchors Excitatory & Astrocyte;\n"
          "Inhibitory anchor omitted — ~0 significant genes). Shown for V1 and V2.\n\n"
          "Statistics caveat: cell-level Wilcoxon DE — adj p-values inflated by\n"
          "pseudoreplication; treat direction and log2FC magnitude as primary.",
          cover_path)
    div1 = tmp / 'd1.pdf'; divider('SET 1 — EtOH-significant genes', div1)
    div2 = tmp / 'd2.pdf'; divider('SET 2 — EtOH-specific genes\n(unchanged by MCT1i)', div2)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    n1 = page_counts['set1']
    out.insert_pdf(fitz.open(str(div1)))
    out.insert_pdf(src, from_page=0, to_page=n1 - 1)
    out.insert_pdf(fitz.open(str(div2)))
    out.insert_pdf(src, from_page=n1, to_page=n1 + page_counts['set2'] - 1)
    out_pdf = args.outdir / f'Xenium_Heatmaps_crossCellType_4contrasts_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, out.page_count if False else 2 + n1 + page_counts['set2'] + 1)

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
