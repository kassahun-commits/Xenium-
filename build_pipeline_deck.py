#!/usr/bin/env python3
"""
Build a GENERAL teaching PowerPoint on running a 10x Xenium analysis with Claude,
plus the statistical caveats to watch for. Broad audience, large fonts.

Sections (confirmed with user 2026-06-09):
  1. The workflow, from the CSV/output files to results (general, large fonts).
  2-3. QC parameters to have Claude check: what each means + recommended threshold;
       plus how to filter by cluster in Xenium Explorer (light steps).
  4. Wilcoxon vs pseudobulk: why/when each, + the "null pseudobulk != no effect"
     diagnostic (effects go the same direction, just underpowered).
  5. The Wilcoxon-vs-pseudobulk significant-gene bar chart.

Regenerates two figures from the V1 whole-cell DE result CSVs, then assembles the
deck. No hardcoded paths (CLI args).
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---- palette -------------------------------------------------------------
NAVY = RGBColor(0x1F, 0x2D, 0x3D)
BLUE = RGBColor(0x1F, 0x77, 0xB4)
RED = RGBColor(0xD6, 0x27, 0x28)
GRAY = RGBColor(0x55, 0x55, 0x55)
LGRAY = RGBColor(0x88, 0x88, 0x88)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GREEN = RGBColor(0x2C, 0xA0, 0x2C)
AMBER = RGBColor(0xE8, 0x8A, 0x00)
PURPLE = RGBColor(0x6B, 0x4E, 0x9E)

# V1 whole-cell contrasts (test, reference) -> pseudobulk comparison string
V1_CONTRASTS = [
    ('EtOH_veh', 'H2O_veh'),
    ('ChronicEtOH', 'H2O_veh'),
    ('H2O_MCT1i', 'H2O_veh'),
    ('EtOH_MCT1i', 'H2O_veh'),
    ('EtOH_MCT1i', 'EtOH_veh'),
    ('MAT2A_OE', 'MAT2A_CM'),
    ('APP', 'H2O_veh'),
]
PB_PADJ, PB_LFC = 0.05, 0.5
WX_PADJ, WX_LFC = 1e-3, 0.5
DIAG_CT, DIAG_TEST, DIAG_REF = 'Neuron', 'EtOH_veh', 'H2O_veh'


# ---- figure helpers ------------------------------------------------------
def _wx_sig(wx):
    return (wx['padj'] < WX_PADJ) & (wx['logfc'].abs() > WX_LFC)


def fig_sig_bar(pb_csv, wx_csv, out_png):
    """Grouped bar: # significant genes per V1 contrast, Wilcoxon vs pseudobulk,
    one panel per cell type."""
    pb = pd.read_csv(pb_csv)
    wx = pd.read_csv(wx_csv)
    wx = wx.copy()
    wx['sig'] = _wx_sig(wx)

    cts = ['Neuron', 'Astrocyte']
    labels = [f'{t}\nvs {r}' for t, r in V1_CONTRASTS]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0))
    for ax, ct in zip(axes, cts):
        pb_vals, wx_vals = [], []
        for t, r in V1_CONTRASTS:
            comp = f'{t}_vs_{r}'
            pb_vals.append(int(pb[(pb['celltype'] == ct) &
                                  (pb['comparison'] == comp) &
                                  (pb['significant'])].shape[0]))
            wx_vals.append(int(wx[(wx['subtype'] == ct) & (wx['test'] == t) &
                                  (wx['reference'] == r) & (wx['sig'])].shape[0]))
        x = np.arange(len(labels))
        w = 0.4
        b1 = ax.bar(x - w / 2, wx_vals, w, color='#d62728',
                    label='Cell-level Wilcoxon (padj<1e-3)')
        b2 = ax.bar(x + w / 2, pb_vals, w, color='#1f77b4',
                    label='Pseudobulk DESeq2 (padj<0.05)')
        ax.set_title(ct, fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_ylabel('# significant genes (|log2FC|>0.5)', fontsize=10)
        ax.bar_label(b1, fontsize=7.5, padding=1)
        ax.bar_label(b2, fontsize=7.5, padding=1)
        for sp in ['top', 'right']:
            ax.spines[sp].set_visible(False)
        if ct == 'Neuron':
            ax.legend(fontsize=9, frameon=False, loc='upper left')
    fig.suptitle('Same cells, two tests — Wilcoxon inflates the hit list',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(out_png, dpi=170, bbox_inches='tight')
    plt.close(fig)


def fig_diagnostic(pb_csv, wx_csv, out_png):
    """Two panels for the EtOH_veh vs H2O_veh neuron contrast:
    (L) pseudobulk volcano showing big-but-underpowered effects below the sig line;
    (R) the diagnostic metrics that say 'underpowered, not null'."""
    pb = pd.read_csv(pb_csv)
    wx = pd.read_csv(wx_csv)
    comp = f'{DIAG_TEST}_vs_{DIAG_REF}'
    d = pb[(pb['celltype'] == DIAG_CT) & (pb['comparison'] == comp)].copy()
    d = d.dropna(subset=['log2FoldChange', 'padj'])

    # metrics
    n_lfc05 = int((d['log2FoldChange'].abs() > 0.5).sum())
    n_lfc1 = int((d['log2FoldChange'].abs() > 1.0).sum())
    min_padj = float(d['padj'].min())
    w = wx[(wx['subtype'] == DIAG_CT) & (wx['test'] == DIAG_TEST) &
           (wx['reference'] == DIAG_REF)].copy()
    wsig = w[_wx_sig(w)]
    m = wsig.merge(d[['gene', 'log2FoldChange']], on='gene')
    conc = (np.sign(m['logfc']) == np.sign(m['log2FoldChange'])).mean() * 100 \
        if len(m) else float('nan')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4),
                             gridspec_kw={'width_ratios': [1.25, 1]})
    # ---- L: volcano ----
    ax = axes[0]
    x = d['log2FoldChange'].values
    y = -np.log10(np.clip(d['padj'].values, 1e-300, 1))
    sig = (d['padj'] < PB_PADJ) & (d['log2FoldChange'].abs() > PB_LFC)
    ax.scatter(x[~sig.values], y[~sig.values], s=7, c='#cccccc', alpha=0.55,
               edgecolor='none', rasterized=True)
    ax.scatter(x[sig.values], y[sig.values], s=12, c='#1f77b4',
               edgecolor='none')
    ax.axhline(-np.log10(PB_PADJ), ls='--', color='k', lw=0.7, alpha=0.5)
    ax.text(ax.get_xlim()[1], -np.log10(PB_PADJ), ' padj=0.05',
            fontsize=8, va='bottom', ha='right', color='k')
    for v in (-PB_LFC, PB_LFC):
        ax.axvline(v, ls='--', color='k', lw=0.7, alpha=0.5)
    # label the largest-effect near-miss genes
    cand = d[d['log2FoldChange'].abs() > 1.0].nsmallest(8, 'padj')
    for _, rr in cand.iterrows():
        ax.annotate(rr['gene'],
                    (rr['log2FoldChange'], -np.log10(max(rr['padj'], 1e-300))),
                    fontsize=8, ha='center', va='bottom', color='#444')
    ax.set_title(f'{DIAG_CT}: {DIAG_TEST} vs {DIAG_REF}\n'
                 f'pseudobulk — 0 pass, but effects are real-sized',
                 fontsize=11)
    ax.set_xlabel('log2 fold change', fontsize=10)
    ax.set_ylabel('-log10 adjusted p', fontsize=10)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)

    # ---- R: metrics ----
    ax = axes[1]
    ax.axis('off')
    lines = [
        ('Is the null real, or underpowered?', '#1F2D3D', 13, True),
        ('', '#1F2D3D', 6, False),
        (f'• {n_lfc05:,} genes with |log2FC| > 0.5', '#333333', 12, False),
        (f'• {n_lfc1:,} genes with |log2FC| > 1.0  (large!)', '#333333', 12, False),
        (f'• best gene only reaches padj = {min_padj:.2f}', '#333333', 12, False),
        ('   (not ~1.0 — knocking on the door)', '#777777', 10, False),
        (f'• {conc:.0f}% of Wilcoxon hits go the SAME', '#1f77b4', 12, True),
        ('   direction in pseudobulk', '#1f77b4', 12, True),
        ('', '#1F2D3D', 8, False),
        ('Verdict: effects are present and consistent,', '#2ca02c', 12, True),
        ('but n=2 vs 2 is too few to certify them.', '#2ca02c', 12, True),
        ('→ A null here means ADD Ns, not "no effect".', '#d62728', 12, True),
    ]
    yy = 0.97
    for text, col, sz, bold in lines:
        if text:
            ax.text(0.0, yy, text, fontsize=sz, color=col,
                    fontweight='bold' if bold else 'normal',
                    transform=ax.transAxes, va='top')
        yy -= 0.075 if text else 0.04
    fig.suptitle('Reading a "no significant DE" pseudobulk result',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(out_png, dpi=170, bbox_inches='tight')
    plt.close(fig)


# ---- pptx helpers --------------------------------------------------------
def add_blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def txt(slide, l, t, w, h, text, size=18, bold=False, color=NAVY,
        align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.color.rgb = color
    r.font.name = 'Calibri'
    return tb


def bullets(slide, l, t, w, h, items, size=18, color=NAVY, space=10):
    tb = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, (lvl, text, *style) in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = lvl
        p.space_after = Pt(space)
        r = p.add_run()
        r.text = ('• ' if lvl == 0 else '– ') + text
        r.font.size = Pt(size - lvl * 2)
        r.font.color.rgb = style[0] if style else color
        r.font.name = 'Calibri'
        if len(style) > 1 and style[1]:
            r.font.bold = True
    return tb


def band(slide, prs, color=NAVY, h=1.05):
    s = slide.shapes.add_shape(1, Inches(0), Inches(0), prs.slide_width,
                               Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = color
    s.line.fill.background()
    s.shadow.inherit = False
    return s


def chip(slide, l, t, w, h, text, fill, tcolor=WHITE, size=13, bold=True):
    s = slide.shapes.add_shape(5, Inches(l), Inches(t), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    s.line.color.rgb = WHITE
    s.line.width = Pt(1)
    s.shadow.inherit = False
    tf = s.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.color.rgb = tcolor
    r.font.name = 'Calibri'
    return s


def arrow(slide, l, t, w, h=0.3, color=LGRAY):
    s = slide.shapes.add_shape(13, Inches(l), Inches(t), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = color
    s.line.fill.background()
    s.shadow.inherit = False
    return s


def slide_header(slide, prs, title, sub=None):
    band(slide, prs)
    txt(slide, 0.5, 0.16, 12.3, 0.62, title, size=28, bold=True, color=WHITE)
    if sub:
        txt(slide, 0.5, 0.74, 12.3, 0.32, sub, size=13,
            color=RGBColor(0xCF, 0xDA, 0xE6))


def style_table(table, header_size=13, body_size=12, bold_col0=True):
    for j, cell in enumerate(table.rows[0].cells):
        cell.fill.solid()
        cell.fill.fore_color.rgb = NAVY
        for p in cell.text_frame.paragraphs:
            for r in p.runs:
                r.font.color.rgb = WHITE
                r.font.bold = True
                r.font.size = Pt(header_size)
                r.font.name = 'Calibri'
    for i in range(1, len(table.rows)):
        for j in range(len(table.columns)):
            c = table.cell(i, j)
            for p in c.text_frame.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(body_size)
                    r.font.color.rgb = NAVY if j == 0 else GRAY
                    r.font.name = 'Calibri'
                    if j == 0 and bold_col0:
                        r.font.bold = True


# ---- build ---------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pb-csv', required=True,
                    help='WholeCell_V1_PseudobulkDE_results CSV')
    ap.add_argument('--wx-csv', required=True,
                    help='PI_Master_NvsA_V1_DE Wilcoxon CSV')
    ap.add_argument('--out', required=True, help='output .pptx path')
    ap.add_argument('--date', default='2026-06-09')
    args = ap.parse_args()

    out = Path(args.out)
    D = args.date
    assets = out.parent / 'pipeline_deck_assets'
    assets.mkdir(exist_ok=True)

    sig_png = assets / 'sig_bar_v1.png'
    diag_png = assets / 'diagnostic_EtOHveh_neuron.png'
    print('[fig] significance bar chart...')
    fig_sig_bar(args.pb_csv, args.wx_csv, sig_png)
    print('[fig] underpowered-vs-null diagnostic...')
    fig_diagnostic(args.pb_csv, args.wx_csv, diag_png)

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    # ---- Slide 1: title ----
    s = add_blank(prs)
    bg = s.shapes.add_shape(1, 0, 0, prs.slide_width, prs.slide_height)
    bg.fill.solid()
    bg.fill.fore_color.rgb = NAVY
    bg.line.fill.background()
    bg.shadow.inherit = False
    txt(s, 0.8, 2.1, 11.7, 1.4, 'Spatial transcriptomics with Claude',
        size=42, bold=True, color=WHITE)
    txt(s, 0.8, 3.5, 11.7, 0.9,
        'A practical 10x Xenium pipeline — and the caveats to watch for',
        size=22, color=RGBColor(0xCF, 0xDA, 0xE6))
    txt(s, 0.8, 6.45, 11.7, 0.5,
        f'{D}   ·   general overview', size=14, color=LGRAY)

    # ---- Slide 2: workflow (general, LARGE) ----
    s = add_blank(prs)
    slide_header(s, prs, 'The workflow: from your files to results')
    steps = [
        ('1\nHand over\nthe files', 'Xenium output folder\n(count matrix +\nmolecule positions)\n+ ROI CSVs from\nExplorer', PURPLE),
        ('2\nBuild the\nmatrix', 'cells × genes table;\nassign each cell to\nan ROI (point-in-\npolygon)', BLUE),
        ('3\nQC\nfilter', 'drop low-quality\ntranscripts, empty\ncells, rare genes', BLUE),
        ('4\nNormalize\n+ log', 'put cells on a\ncomparable scale', BLUE),
        ('5\nCell\ntyping', 'score marker genes\n→ Neuron, Astrocyte,\nOligo, Microglia', BLUE),
        ('6\nDifferential\nexpression', 'compare groups;\nvolcano plots +\nresult tables', GREEN),
    ]
    x = 0.55
    bw = 1.95
    gap = 0.12
    y = 1.9
    for i, (head, body, col) in enumerate(steps):
        chip(s, x, y, bw, 1.55, head, col, size=15)
        txt(s, x, y + 1.62, bw, 1.7, body, size=12.5, color=GRAY)
        if i < len(steps) - 1:
            arrow(s, x + bw + 0.002, y + 0.62, gap + 0.02)
        x += bw + gap
    txt(s, 0.55, 5.9, 12.2, 1.3,
        'Key idea: Xenium does not hand you "the data" as a table of genes per cell '
        'type. It gives molecule positions and a raw count matrix. Everything else — '
        'ROIs, cell types, expression per group — is computed downstream. You give the '
        'files and the questions; the pipeline builds the rest.',
        size=16, color=NAVY)

    # ---- Slide 3: downstream QC (analyst's job) ----
    s = add_blank(prs)
    slide_header(s, prs, 'Downstream QC — the steps left to you',
                 'The instrument already filtered transcript quality (QV>=20). '
                 'These are what you (or Claude) still need to do.')
    rows = [
        ('Downstream QC step', 'Brain-tissue starting point', 'Why / how to set it'),
        ('Low-count cells', 'counts/cell >= ~10',
         'Drop empty / debris segmentations. Keep permissive — small glia '
         '(microglia) are genuinely low-count.'),
        ('Low-gene cells', 'genes/cell >= ~5',
         'Removes near-empty cells with almost nothing detected.'),
        ('Rare genes', 'gene in >= 5 cells',
         'Drops genes seen in almost no cells (unstable, mostly noise).'),
        ('Doublets / over-segmentation', 'flag very high counts or large area',
         'Dense brain + neuropil merges neighbouring cells into one segment.'),
        ('Cell-area outliers', 'drop implausibly tiny / huge cells',
         'Catches segmentation errors common in tightly packed tissue.'),
        ('Negative-control rate', 'check probe + codeword rate (per sample)',
         'Estimates ambient false-positive level. No fixed cutoff — compare across samples.'),
    ]
    gt = s.shapes.add_table(len(rows), 3, Inches(0.5), Inches(1.5),
                            Inches(12.3), Inches(4.1)).table
    gt.columns[0].width = Inches(3.4)
    gt.columns[1].width = Inches(3.1)
    gt.columns[2].width = Inches(5.8)
    for i in range(len(rows)):
        for j in range(3):
            gt.cell(i, j).text = rows[i][j]
    style_table(gt, header_size=14, body_size=12.5)
    for i in range(1, len(rows)):
        for p in gt.cell(i, 1).text_frame.paragraphs:
            for r in p.runs:
                r.font.bold = True
                r.font.color.rgb = BLUE
    txt(s, 0.5, 5.8, 12.3, 1.4,
        'None of these are fixed 10x numbers — 10x says plot each distribution and set '
        'cutoffs from the data (e.g. 3-5x MAD from the median), starting permissive and '
        'keeping them identical across samples. Do it in Python (prompt Claude / scanpy) '
        'or yourself in R (Seurat / Bioconductor).',
        size=14, color=NAVY)

    # ---- Slide 4: filter by cluster ----
    s = add_blank(prs)
    slide_header(s, prs, 'Optional: pull just one cluster / cell population',
                 'When you only want the cells in cluster X')
    chip(s, 0.55, 1.7, 5.9, 0.55, 'In Xenium Explorer', PURPLE, size=15)
    bullets(s, 0.6, 2.45, 6.0, 3.6, [
        (0, 'Open the Cell groups panel (right side).', NAVY),
        (0, 'Click the cluster you want → "Select cells in group".', NAVY),
        (0, 'Open the Selection panel → download as CSV.', NAVY),
        (0, 'You get cell IDs + coordinates for that cluster only.', GRAY),
    ], size=16)
    chip(s, 6.85, 1.7, 5.9, 0.55, 'Or filter afterward (robust)', BLUE, size=15)
    bullets(s, 6.9, 2.45, 5.9, 3.6, [
        (0, 'Export the full clustering CSV (all cells + cluster column).', NAVY),
        (0, 'Filter to your cluster in pandas:', NAVY),
        (1, 'df[df["Cluster"] == 1]', GRAY),
        (0, 'Works regardless of Explorer version.', GRAY),
    ], size=16)
    txt(s, 0.55, 6.0, 12.2, 1.2,
        'Heads-up: an Explorer cluster/selection CSV contains cell IDs + coordinates, '
        'NOT expression counts. To analyze expression for those cells, hand the cell-ID '
        'list to Claude and subset the count matrix on them.',
        size=15, color=AMBER, bold=True)

    # ---- Slide 5: Wilcoxon vs pseudobulk (why each) ----
    s = add_blank(prs)
    slide_header(s, prs, 'Two ways to test differential expression',
                 'Pick the test that matches your replicate structure')
    chip(s, 0.55, 1.55, 5.95, 0.55, 'Cell-level Wilcoxon', RED, size=15)
    bullets(s, 0.6, 2.25, 6.0, 4.6, [
        (0, 'Treats every CELL as a sample.', NAVY, True),
        (0, 'Fast; good for marker discovery and exploration.', NAVY),
        (0, 'Reasonable when you truly have n=1 (e.g. one punch).', NAVY),
        (0, 'BUT with 2–3 mice it is pseudoreplication: p-values track '
            'cell COUNT, not biology → inflated hit lists.', RED, True),
        (0, 'Reviewers discount these p-values.', GRAY),
    ], size=15)
    chip(s, 6.85, 1.55, 5.95, 0.55, 'Pseudobulk DESeq2', BLUE, size=15)
    bullets(s, 6.9, 2.25, 5.95, 4.6, [
        (0, 'Sums counts per MOUSE, then tests across animals.', NAVY, True),
        (0, 'Mouse = the real experimental unit.', NAVY),
        (0, 'Rigorous and reviewer-preferred when n≥2–3.', BLUE, True),
        (0, 'Conservative: needs the effect to reproduce across animals.', NAVY),
        (0, 'Use this as the PRIMARY analysis.', BLUE, True),
    ], size=15)

    # ---- Slide 6: the diagnostic (null != no effect) ----
    s = add_blank(prs)
    slide_header(s, prs, 'A "no significant DE" result is not always "no effect"',
                 'Worked example: EtOH_veh vs H2O_veh, neurons (n=2 vs 2)')
    s.shapes.add_picture(str(diag_png), Inches(0.7), Inches(1.5),
                         width=Inches(12.0))
    txt(s, 0.55, 6.55, 12.2, 0.85,
        'Hundreds of large fold changes, all pointing the same way as the Wilcoxon '
        'result, capped just under significance by sample size. That is an '
        'underpowered effect — the fix is more animals, not a looser test.',
        size=14, color=NAVY)

    # ---- Slide 7: the bar chart ----
    s = add_blank(prs)
    slide_header(s, prs, 'Same data, two tests — the size of the inflation',
                 'Significant genes per contrast (|log2FC|>0.5), V1 whole-cell')
    s.shapes.add_picture(str(sig_png), Inches(0.6), Inches(1.45),
                         width=Inches(12.1))
    txt(s, 0.55, 6.35, 12.2, 1.0,
        'Wilcoxon (red) returns far more "hits" than pseudobulk (blue) — even though '
        'pseudobulk uses the looser p-threshold here. The gap is pseudoreplication. '
        'The contrasts that survive pseudobulk (e.g. the MAT2A and EtOH_MCT1i comparisons) '
        'are the ones with the strongest, most reproducible effects.',
        size=13.5, color=NAVY)

    # ---- Slide 8: takeaways ----
    s = add_blank(prs)
    slide_header(s, prs, 'Takeaways')
    bullets(s, 0.6, 1.55, 12.2, 4.8, [
        (0, 'You provide the files + the question; the pipeline builds matrix → ROIs '
            '→ QC → cell types → DE.', NAVY, True),
        (0, 'The machine filters transcript quality (QV>=20) for you; you still must QC '
            'cells and genes — set cutoffs from the data and keep them identical across samples.', NAVY),
        (0, 'Match the DE test to your design: pseudobulk when you have biological '
            'replicates; Wilcoxon only for exploration / true n=1.', BLUE, True),
        (0, 'A null pseudobulk result with large, consistent fold changes = underpowered '
            '→ add Ns. Fold changes hugging zero = likely truly no effect.', NAVY),
        (0, 'Watch for confounds: if a treatment lives on only one slide/batch, its "hits" '
            'may be batch, not biology.', AMBER, True),
        (0, 'Every figure ships with editable text + a source-data CSV; scripts are version-'
            'controlled.', GRAY),
    ], size=17)

    prs.save(str(out))
    print(f'[saved] {out}')
    print(f'[assets] {sig_png.name}, {diag_png.name}')


if __name__ == '__main__':
    main()
