#!/usr/bin/env python3
"""
Build a single summary PDF for the PI containing:
  - Cover page (study design, slides, sample counts, methodology, caveats)
  - Figure 1: "Everything" curated dot plot (118 genes)
  - Figure 2: Focused dot plot (filtered + force-included, 65 genes)
  - Figure 3: Box-and-whisker plot for the same 65 focused genes
  - (Bonus) Figure 4: Per-mouse mean bar plot for the same 65 focused genes
"""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz  # PyMuPDF

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()

OUT = BASE / f'PI_Summary_XeniumAmygdala_SlidesAB_{DATE}.pdf'

# Source figures (already generated)
EVERYTHING_DOT = BASE / 'SlidesAB_CuratedDotPlot_NeuronAstrocyte_2026-05-24.pdf'
FOCUSED_DOT    = BASE / 'SlidesAB_FocusedDotPlot_v2_NeuronAstrocyte_2026-05-24.pdf'
FOCUSED_BOX    = BASE / 'SlidesAB_FocusedBoxPlots_NeuronAstrocyte_2026-05-24.pdf'
FOCUSED_BAR    = BASE / 'SlidesAB_FocusedBarPlot_v2_NeuronAstrocyte_2026-05-24.pdf'

# ----- Build cover page -----
cover = fitz.open()
page = cover.new_page(width=612, height=792)  # US Letter

# Heading
H = fitz.Rect(54, 60, 558, 105)
page.insert_textbox(H,
                    "Xenium spatial-transcriptomic profiling of mouse hippocampus\n"
                    "Preliminary cell-type-resolved gene expression — Summary for PI",
                    fontsize=15, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 108, 558, 124),
                    f"Naomi Kassahun · {DATE} · MEWS Lab",
                    fontsize=10, fontname='Helvetica-Oblique')

# Section: Study design
y = 145
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Study design", fontsize=12, fontname='Helvetica-Bold')
y += 22
design_text = (
    "Two Xenium slides (Cassette A, Cassette B) of mouse hippocampus tissue punches.\n"
    "Panel: Mews Single Cell-1 (98 custom genes) on Xenium Mouse 5K Pan-Tissue base "
    "panel (~5,000 genes total).\n"
    "Tissue punches assigned to experimental groups via ROIs drawn in Xenium Explorer; "
    "cells assigned to ROIs by point-in-polygon on cell centroids."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 60),
                    design_text, fontsize=9.5, fontname='Helvetica')

# Section: Sample composition
y = 215
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Samples after combining Slides A + B (n = mice per group)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
samples_text = (
    "    H2O + veh          n = 2 mice            EtOH + MCT1i    n = 3 mice\n"
    "    H2O + MCT1i      n = 2 mice            Chronic EtOH      n = 3 mice\n"
    "    EtOH + veh          n = 2 mice            MAT2A_CM          n = 3 mice\n"
    "                                                       MAT2A_OE           n = 3 mice\n"
    "Total: 167,677 cells passed QC across 18 samples in 7 groups.\n"
    "APP punches not included (not drawn in Xenium Explorer for this round)."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 100),
                    samples_text, fontsize=9.5, fontname='Helvetica')

# Section: Cell typing
y = 320
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Cell-type assignment", fontsize=12, fontname='Helvetica-Bold')
y += 22
ct_text = (
    "Marker-module scoring (scanpy.tl.score_genes) against four lineages:\n"
    "  Neuron:  Rbfox3, Snap25, Syn1, Syt1, Stmn2, Map2, Tubb3\n"
    "  Astrocyte:  Gfap, Aqp4, Slc1a3, Aldh1l1, S100b, Aldoc, Gja1\n"
    "  Oligodendrocyte:  Mog, Olig1/2, Sox10\n"
    "  Microglia:  Cx3cr1, Tmem119, Csf1r, Aif1, Trem2\n"
    "Each cell assigned to highest-scoring lineage (score > 0); otherwise \"Unclassified\".\n"
    "Marker validation: every assigned cell type expresses its own markers highly\n"
    "and other lineages' markers weakly — confirmed in a prior validation plot."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 130),
                    ct_text, fontsize=9.5, fontname='Helvetica')

# Section: Figures in this PDF
y = 460
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Figures in this PDF", fontsize=12, fontname='Helvetica-Bold')
y += 22
figs_text = (
    "Figure 1 (next page).  All curated genes (118 genes from 5K panel), dot plot.\n"
    "  Categories: cell-type identity, alcohol/acetaldehyde metab, acetate/Ac-CoA, lactate/MCT,\n"
    "  glucose transport, glycolysis/TCA, SAM/MAT2A, glutathione/oxidative stress, stress/UPR,\n"
    "  activity/IEGs, glutamate signaling, GABA signaling, endocannabinoid, neuroinflammation,\n"
    "  glutamate transporters, HATs (lysine acetyltransferases), HDACs / Sirtuins.\n"
    "  Color = mean log-norm expression. Dot size = fraction of cells expressing.\n\n"
    "Figure 2.  Focused dot plot (65 genes). Same data, filtered to genes that visually vary\n"
    "  either between Neurons and Astrocytes (|Δ| ≥ 0.5) or between EtOH+veh and EtOH+MCT1i\n"
    "  (|Δ| ≥ 0.3); plus key axis genes force-included (MCTs, GLUTs, ACSS1/2, MAT2A axis,\n"
    "  p300/CBP/KAT2A/KAT5, HDAC1/2/3, SIRT1, SLC1A1/3).\n\n"
    "Figure 3.  Box-and-whisker plot for the same 65 focused genes.\n"
    "  One panel per gene; 14 boxes per panel (7 condition groups × Neuron/Astrocyte).\n"
    "  Boxes: cell-level IQR + median; whiskers: 1.5×IQR; outliers hidden; black dots = per-mouse means.\n\n"
    "Figure 4 (bonus).  Per-mouse mean bar plot of the same 65 genes — clearer view of\n"
    "  between-mouse variability for sparse data (n=2-3 biological replicates per group)."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 220),
                    figs_text, fontsize=9, fontname='Helvetica')

# Caveats
y = 690
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Caveats", fontsize=12, fontname='Helvetica-Bold')
y += 22
caveat_text = (
    "Preliminary. n = 2 or 3 mice per group — figures show qualitative/trend-level patterns.\n"
    "Xenium spatial transcript assignment can produce some cross-cell-type contamination\n"
    "(transcript bleed) — e.g., neuronal markers occasionally appearing weakly in astrocyte\n"
    "panels. Differential expression p-values (when reported separately) are exploratory\n"
    "due to cell-level pseudoreplication; pseudobulk would be required for publication-grade\n"
    "statistics."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    caveat_text, fontsize=8.5, fontname='Helvetica-Oblique')

# Save cover
COVER_PATH = BASE / '_PI_Summary_cover_tmp.pdf'
cover.save(str(COVER_PATH))
cover.close()
print(f'Cover -> {COVER_PATH}')

# ----- Merge cover + figure PDFs -----
out = fitz.open()
out.insert_pdf(fitz.open(str(COVER_PATH)))
for label, path in [('Figure 1 — Everything dot plot', EVERYTHING_DOT),
                    ('Figure 2 — Focused dot plot', FOCUSED_DOT),
                    ('Figure 3 — Focused box-and-whisker', FOCUSED_BOX),
                    ('Figure 4 — Focused bar plot (bonus)', FOCUSED_BAR)]:
    if not path.exists():
        print(f'MISSING: {path}'); continue
    out.insert_pdf(fitz.open(str(path)))
    print(f'  added {label}: {path.name}')

out.save(str(OUT))
out.close()
print(f'\nWrote {OUT}')

# Clean up tmp cover
os.remove(COVER_PATH)
