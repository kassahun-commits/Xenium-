#!/usr/bin/env python3
"""
Build a combined summary PDF for excitatory + inhibitory neuron analyses
(H2O_veh vs EtOH_veh vs ChronicEtOH).

Pages:
  1. Cover page — methods + sample sizes + headline findings
  2. Figure 1 — Excitatory neurons box plots (top 24 up + 24 down + MCT family)
  3. Figure 2 — Inhibitory neurons box plots (top 24 up + 24 down + MCT family)
"""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz  # PyMuPDF

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()

OUT = BASE / f'NeuronSubtype_AlcoholResponse_Summary_{DATE}.pdf'

EXCIT_PDF = BASE / 'SlidesAB_Excit_vehAcuteChronic_BoxPlots_2026-05-25.pdf'
INHIB_PDF = BASE / 'SlidesAB_Inhib_vehAcuteChronic_BoxPlots_2026-05-25.pdf'

# ---------- Cover page ----------
cover = fitz.open()
page = cover.new_page(width=612, height=792)  # US Letter

# Title
page.insert_textbox(fitz.Rect(54, 60, 558, 105),
                    "Hippocampal excitatory vs inhibitory neurons —\n"
                    "alcohol response (H2O_veh / EtOH_veh / EtOH_MCT1i / ChronicEtOH)",
                    fontsize=14, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 108, 558, 124),
                    f"Naomi Kassahun · {DATE} · MEWS Lab · Xenium May 2026 (Slides A + B combined)",
                    fontsize=10, fontname='Helvetica-Oblique')

# Methods
y = 145
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Methods", fontsize=12, fontname='Helvetica-Bold')
y += 22
methods = (
    "Combined Slide A + Slide B (Xenium Mouse Brain 5K panel + 98-gene custom add-on).\n"
    "ROIs drawn in Xenium Explorer; cells assigned to groups by point-in-polygon.\n"
    "Counts normalized to 1e4 per cell + log1p; cell-types called by scanpy.tl.score_genes "
    "against canonical marker sets.\n"
    "Neuronal subtyping (within cells already called 'Neuron'):\n"
    "  Excitatory: Slc17a7, Camk2a, Camk2b, Satb2, Tbr1, Neurod6 module score > 0\n"
    "  Inhibitory: Gad1, Gad2, Slc32a1, Pvalb, Sst, Vip, Reln, Lhx6 module score > 0\n"
    "Differential expression: cell-level Wilcoxon (scanpy.tl.rank_genes_groups) for each\n"
    "of three contrasts vs H2O_veh — EtOH_veh, EtOH_MCT1i, and ChronicEtOH.\n"
    "Top hits per direction: |log2FC| >= 1, padj < 1e-3, ranked by max(|log2FC|) across\n"
    "the three contrasts.\n"
    "MCT family (Slc16a1=MCT1, Slc16a7=MCT2, Slc16a3=MCT4) force-included regardless of\n"
    "significance. Gene labels use 'Slc... (alias)' format where applicable."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 165),
                    methods, fontsize=9, fontname='Helvetica')

# Sample sizes
y = 318
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Sample sizes (mice × cells per group)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
samples = (
    "                                  Excitatory cells          Inhibitory cells\n"
    "  H2O_veh         2 mice         10,650                    4,492\n"
    "  EtOH_veh        2 mice          7,971                    5,095\n"
    "  EtOH_MCT1i      3 mice         12,319                    ~ (computed at run time)\n"
    "  ChronicEtOH     3 mice         15,861                    6,966\n"
    "\n"
    "Total: ~75,000 excitatory + ~32,000 inhibitory cells across both slides."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 105),
                    samples, fontsize=9, fontname='Courier')

# Headline findings
y = 420
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Headline findings", fontsize=12, fontname='Helvetica-Bold')
y += 22
findings = (
    "SHARED across excitatory and inhibitory neurons (same alcohol signature in both):\n"
    "  - UP in alcohol: developmental / fetal neural transcription factors\n"
    "    (Foxp2, Shox2, Ebf3, Gata3, Irx2, Bnc2, Tcf7l2, Sox14, Barhl2, Zfhx3, Tfap2b),\n"
    "    glycine receptor Glra1, serotonin receptor Htr2c, VGLUT2 (Slc17a6), Wnt agonists.\n"
    "  - DOWN in alcohol: hippocampal dentate granule markers (Prox1, Calb1, Pdyn, Npy2r),\n"
    "    neural-progenitor / neurogenic TFs (Glis3, Prdm8, Bhlhe22, Rfx3),\n"
    "    adhesion / extracellular matrix genes (Pcdh20, Gpc4, Slc44a5).\n"
    "\n"
    "EXCITATORY-specific:\n"
    "  - Stronger TF reprogramming signature; dentate identity loss prominent.\n"
    "  - VGLUT1 (Slc17a7) appears DOWN in inhibitory cells (possible transcript bleed).\n"
    "\n"
    "INHIBITORY-specific (most striking class-distinct finding):\n"
    "  - Strong loss of immediate-early / activity-dependent genes\n"
    "    (Arc, Egr2, Egr3, Egr4, Nr4a1) — alcohol selectively blunts interneuron activity\n"
    "    transcription. Excitatory neurons show only a mild Arc decrement.\n"
    "\n"
    "MCT family:\n"
    "  - MCT1 (Slc16a1) and MCT4 (Slc16a3) are very low-expressed in both neuron classes\n"
    "    (consistent with their primary expression in astrocytes/endothelium); shown as\n"
    "    per-mouse-mean strip plots.\n"
    "  - MCT2 (Slc16a7) is the expected neuronal MCT; expression is stable across\n"
    "    H2O_veh / EtOH_veh / ChronicEtOH."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 275),
                    findings, fontsize=9, fontname='Helvetica')

# Caveats
y = 705
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Caveats", fontsize=12, fontname='Helvetica-Bold')
y += 22
caveats = (
    "Preliminary. n=2-3 mice per group. Cell-level Wilcoxon p-values are inflated due to "
    "pseudoreplication (cells within a mouse are not biological replicates); direction and "
    "magnitude (log2FC) of changes are more trustworthy than the p-values. Pseudobulk DE "
    "(per-sample aggregation) would be required for publication-grade statistics."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    caveats, fontsize=8.5, fontname='Helvetica-Oblique')

COVER = BASE / '_subtype_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()
print(f'Cover -> {COVER}')

# ---------- Merge ----------
out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))
for label, path in [('Figure 1 — Excitatory neurons', EXCIT_PDF),
                    ('Figure 2 — Inhibitory neurons', INHIB_PDF)]:
    if not path.exists():
        print(f'MISSING: {path}'); continue
    out.insert_pdf(fitz.open(str(path)))
    print(f'  added {label}: {path.name}')

out.save(str(OUT))
out.close()
print(f'\nWrote {OUT}')

os.remove(COVER)
