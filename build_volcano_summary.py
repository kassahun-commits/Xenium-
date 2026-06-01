#!/usr/bin/env python3
"""Build a cover + 18-volcano summary PDF."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz

import sys
BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()
# Optional CLI: pass lfc threshold + suffix to label outputs
LFC = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
SUFFIX = sys.argv[2] if len(sys.argv) > 2 else ''
OUT = BASE / f'NeuronAstrocyte_VolcanoGrid_Summary_{DATE}{SUFFIX}.pdf'
SRC = BASE / f'SlidesAB_volcanos_AllContrasts_{DATE}{SUFFIX}.pdf'

cover = fitz.open()
page = cover.new_page(width=612, height=792)
page.insert_textbox(fitz.Rect(54, 60, 558, 110),
                    "Hippocampal Excitatory | Inhibitory | Astrocyte —\n"
                    "differential expression volcano plots across 6 contrasts",
                    fontsize=14, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 112, 558, 128),
                    f"Naomi Kassahun · {DATE} · MEWS Lab · Xenium May 2026 (Slides A + B)",
                    fontsize=10, fontname='Helvetica-Oblique')

y = 150
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Layout", fontsize=12, fontname='Helvetica-Bold')
y += 22
layout = (
    "One comparison per page; three volcano panels per page:\n"
    "  Excitatory (blue header)  |  Inhibitory (purple header)  |  Astrocyte (red header)\n"
    "Red dots = up in the test group. Blue dots = up in the reference group.\n"
    "Gray dots = non-significant. Dashed lines mark the significance thresholds.\n"
    "Top ~12 genes per direction are labeled by name.\n"
    "Each panel has a small inset showing: test/reference cell counts, mouse counts,\n"
    "and the number of significant up- and down-regulated genes.\n"
    "\n"
    f"Significance thresholds: |log2 fold change| >= {LFC:g} AND adjusted p < 0.001."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 160),
                    layout, fontsize=9, fontname='Helvetica')

y = 320
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Contrasts (one per page)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
contrasts = (
    "  Page 1.  EtOH_veh        vs H2O_veh        — Acute alcohol\n"
    "  Page 2.  ChronicEtOH     vs H2O_veh        — Chronic alcohol\n"
    "  Page 3.  H2O_MCT1i       vs H2O_veh        — Drug only (does MCT1i alone do anything?)\n"
    "  Page 4.  EtOH_MCT1i      vs H2O_veh        — Alcohol + drug vs control (combined effect)\n"
    "  Page 5.  EtOH_MCT1i      vs EtOH_veh       — MCT1i rescue effect (within alcohol-treated)\n"
    "  Page 6.  MAT2A_OE        vs MAT2A_CM       — MAT2A overexpression\n"
    "                                              (vs catalytic-mutant control, NOT vs H2O_veh)"
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 160),
                    contrasts, fontsize=9, fontname='Courier')

y = 500
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Statistics",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
stats = (
    "Cell-level Wilcoxon rank-sum test (scanpy.tl.rank_genes_groups).\n"
    "Tested independently per cell-type subset per contrast.\n"
    "Adjusted p-values: Benjamini-Hochberg.\n"
    "\n"
    "CAVEAT: cell-level Wilcoxon with n=2-3 mice per group is pseudoreplicated.\n"
    "Adjusted p-values are inflated by the large per-mouse cell counts.\n"
    "Treat the direction and magnitude (log2FC) of changes as reliable;\n"
    "treat the adjusted p-values as ranking statistics rather than as biological\n"
    "evidence of replicability. For publication-grade statistics, run pseudobulk\n"
    "DE (per-sample aggregation, then DESeq2 / limma-voom across samples)."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 200),
                    stats, fontsize=9, fontname='Helvetica')

y = 700
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Cell-type definitions",
                    fontsize=10, fontname='Helvetica-Bold')
y += 16
ct = (
    "Excitatory: Neuron lineage + Slc17a7/Camk2a/b/Satb2/Tbr1/Neurod6 score >0 (~75,000 cells).\n"
    "Inhibitory: Neuron lineage + Gad1/Gad2/Slc32a1/Pvalb/Sst/Vip/Reln/Lhx6 score >0 (~32,000 cells).\n"
    "Astrocyte: top-scoring lineage = Astrocyte (Gfap/Aqp4/Slc1a3/Aldh1l1/...) (~44,000 cells)."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    ct, fontsize=8, fontname='Helvetica-Oblique')

COVER = BASE / '_volcano_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()

out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))
out.insert_pdf(fitz.open(str(SRC)))
out.save(str(OUT))
out.close()
print(f'Wrote {OUT}')
os.remove(COVER)
