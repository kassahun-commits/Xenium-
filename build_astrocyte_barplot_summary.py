#!/usr/bin/env python3
"""Build a combined PDF for the astrocyte bar-plot grid (cover + 7 pages)."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()
OUT = BASE / f'Astrocyte_AllGroups_BarPlot_Summary_{DATE}.pdf'

ASTRO = BASE / f'SlidesAB_Astrocyte_AllGroups_BarPlots_{DATE}.pdf'

# Cover page
cover = fitz.open()
page = cover.new_page(width=612, height=792)
page.insert_textbox(fitz.Rect(54, 60, 558, 110),
                    "Hippocampal astrocytes — per-mouse mean bar plots\n"
                    "across all 7 condition groups",
                    fontsize=14, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 112, 558, 128),
                    f"Naomi Kassahun · {DATE} · MEWS Lab · Xenium May 2026 (Slides A + B)",
                    fontsize=10, fontname='Helvetica-Oblique')

y = 150
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Methods", fontsize=12, fontname='Helvetica-Bold')
y += 22
methods = (
    "Combined Slide A + Slide B (Xenium Mouse 5K + 98-gene Mews custom add-on).\n"
    "ROIs from Xenium Explorer; cells assigned to groups by point-in-polygon.\n"
    "Counts normalized to 1e4 per cell + log1p.\n"
    "Cell-types: scanpy.tl.score_genes vs canonical markers.\n"
    "Astrocyte subset: cells whose top-scoring lineage is the Astrocyte module\n"
    "(Gfap, Aqp4, Slc1a3, Aldh1l1, S100b, Aldoc, Gja1) with positive score.\n"
    "No sub-classification — all astrocytes pooled.\n"
    "\n"
    "Bars = mean of per-mouse means. Error bars = SEM across mice.\n"
    "Black dots = individual mouse means. Gene labels use 'gene (alias)' format\n"
    "where applicable; reporter is labeled 'GFP (reporter)'."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 170),
                    methods, fontsize=9, fontname='Helvetica')

y = 330
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Sample sizes",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
samples = (
    "                              Astrocyte cells\n"
    "  H2O_veh         2 mice\n"
    "  H2O_MCT1i       2 mice\n"
    "  EtOH_veh        2 mice\n"
    "  EtOH_MCT1i      3 mice\n"
    "  ChronicEtOH     3 mice\n"
    "  MAT2A_CM        3 mice\n"
    "  MAT2A_OE        3 mice\n"
    "Total: 44,085 astrocyte cells across 18 samples."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 130),
                    samples, fontsize=9, fontname='Courier')

y = 480
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Gene list (128 genes plotted; 20 per page over 7 pages)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
genes_text = (
    "Astrocyte identity (sanity): Gfap, Aqp4, Slc1a3 (GLAST), Aldh1l1, S100b,\n"
    "Aldoc, Gja1, Slc1a1 (EAAC1), Slc6a11 (GAT-3).\n"
    "Neuronal contrast: Rbfox3 (NeuN), Snap25, Syn1.\n"
    "MCT family: Slc16a1 (MCT1), Slc16a7 (MCT2), Slc16a3 (MCT4).\n"
    "Reporter: GFP.\n"
    "Alcohol-metabolizing enzymes (NEW): Adh1, Adh5, Aldh2, Aldh1a2, Cyp2e1,\n"
    "Cat (catalase).\n"
    "Acetate/Ac-CoA + glycolysis/TCA: Acss1, Acss2, Acly, Eno1, Eno2, Ldha, Ldhb.\n"
    "Oxidative stress/detox: Sod1, Sod2, Gpx1, Gpx4, Nfe2l2, Hmox1.\n"
    "Chromatin / methylation / one-carbon: Ahcy, Brd9, Chd3, Chd8, Chdh, Dpf1,\n"
    "Ino80, Kat7, Kmt5a, Mat1a, Mat2a, Mat2b, Mbd1, Mbd4, Mtap, Mthfd1, Mthfd2,\n"
    "Mtr, Mtrr, Nsd2, Nsd3, Pemt, Phgdh, Prdm2, Prdm8, Prmt2/6/8/9, Psph,\n"
    "Setdb2, Shmt2, Slc25a32, Ss18, Suv39h2, Tyms, Gnmt.\n"
    "IEGs / activity: Atf3, Btg2, Fos, Fosb, Gadd45b, Jun, Junb.\n"
    "Synaptic / RNA-binding: Grik5, Rgs14, Pcp4, Cacng5, Cacnb4, Dlg2, Lrrtm2/4,\n"
    "Pclo, Unc13a, Nptx1, Tagln3, Celf2, Ccnt2, Cdk2ap1, Hmgb3, Hmgn3, Ilf3,\n"
    "Kcnh7, Malat1, Mbnl2, Pnisr, Prpf39, Rbm34, Rsrc2, Slc25a40, Tia1, Zranb2.\n"
    "Other metabolism: Slc5a8, Emb, Bsg, Slc4a4, Crat, Crot, Slc10a2, Acat1, Acat2,\n"
    "Gapdh, Slc13a5, Pdha1, Pdhb, Glul (glutamine synthetase), Abcc3, Abhd14b,\n"
    "Pfkm, Cfhr1, Dgkb, Gng12, Elovl5."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 240),
                    genes_text, fontsize=8.5, fontname='Helvetica')

y = 730
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Caveats", fontsize=10, fontname='Helvetica-Bold')
y += 16
caveats = (
    "Preliminary. n=2-3 mice per group. Bars show per-mouse mean ± SEM; "
    "dots = individual mice. No SEM if n=1."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    caveats, fontsize=8, fontname='Helvetica-Oblique')

COVER = BASE / '_astro_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()

out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))
out.insert_pdf(fitz.open(str(ASTRO)))
out.save(str(OUT))
out.close()
print(f'Wrote {OUT}')
os.remove(COVER)
