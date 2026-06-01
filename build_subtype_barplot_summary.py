#!/usr/bin/env python3
"""Build a combined PDF for excitatory + inhibitory subtype bar-plot grids
(7 pages each + cover page)."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()
OUT = BASE / f'NeuronSubtype_AllGroups_BarPlot_Summary_{DATE}.pdf'

EXCIT = BASE / f'SlidesAB_Excitatory_AllGroups_BarPlots_{DATE}.pdf'
INHIB = BASE / f'SlidesAB_Inhibitory_AllGroups_BarPlots_{DATE}.pdf'

# --- Cover page ---
cover = fitz.open()
page = cover.new_page(width=612, height=792)
page.insert_textbox(fitz.Rect(54, 60, 558, 110),
                    "Hippocampal excitatory & inhibitory neurons —\n"
                    "per-mouse mean bar plots across all 7 condition groups",
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
    "Excitatory subset: Slc17a7, Camk2a/b, Satb2, Tbr1, Neurod6 module score > 0.\n"
    "Inhibitory subset: Gad1, Gad2, Slc32a1, Pvalb, Sst, Vip, Reln, Lhx6 module score > 0.\n"
    "\n"
    "Bars = mean of per-mouse means. Error bars = SEM across mice.\n"
    "Black dots = individual mouse means (so you can see the n=2-3 replicates).\n"
    "Gene labels use 'Slc... (protein alias)' format where applicable; reporter\n"
    "is labeled 'GFP (reporter)'."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 165),
                    methods, fontsize=9, fontname='Helvetica')

y = 325
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Sample sizes",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
samples = (
    "                              Excitatory cells          Inhibitory cells\n"
    "  H2O_veh         2 mice    10,650                    4,492\n"
    "  H2O_MCT1i       2 mice    12,752                    (computed at runtime)\n"
    "  EtOH_veh        2 mice     7,971                    5,095\n"
    "  EtOH_MCT1i      3 mice    12,319                    (computed at runtime)\n"
    "  ChronicEtOH     3 mice    15,861                    6,966\n"
    "  MAT2A_CM        3 mice     ~10,043                  (computed at runtime)\n"
    "  MAT2A_OE        3 mice     ~5,444                   (computed at runtime)"
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 130),
                    samples, fontsize=9, fontname='Courier')

y = 480
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Gene list (122 genes plotted; 20 per page)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
genes_text = (
    "Cell-type identity (for sanity-check rows): Rbfox3 (NeuN), Snap25, Syn1, Map2,\n"
    "Slc17a7 (VGLUT1), Slc17a6 (VGLUT2), Camk2a, Camk2b, Satb2, Tbr1, Neurod6,\n"
    "Gad1, Gad2, Slc32a1 (VGAT), Pvalb, Sst, Vip, Reln, Gfap, Aqp4, Slc1a3 (GLAST),\n"
    "Slc1a1 (EAAC1).\n"
    "MCT family: Slc16a1 (MCT1), Slc16a7 (MCT2), Slc16a3 (MCT4).\n"
    "Reporter: GFP.\n"
    "Chromatin / methylation / one-carbon: Ahcy, Brd9, Chd3, Chd8, Chdh, Dpf1, Ino80,\n"
    "Kat7, Kmt5a, Mat1a, Mat2a, Mat2b, Mbd1, Mbd4, Mtap, Mthfd1, Mthfd2, Mtr, Mtrr,\n"
    "Nsd2, Nsd3, Pemt, Phgdh, Prdm2, Prdm8, Prmt2/6/8/9, Psph, Setdb2, Shmt2,\n"
    "Slc25a32, Ss18, Suv39h2, Tyms, Gnmt.\n"
    "IEGs / activity: Atf3, Btg2, Fos, Fosb, Gadd45b, Jun, Junb.\n"
    "Synaptic / RNA-binding: Grik5, Rgs14, Pcp4, Cacng5, Cacnb4, Dlg2, Lrrtm2/4,\n"
    "Pclo, Unc13a, Nptx1, Tagln3, Celf2, Ccnt2, Cdk2ap1, Hmgb3, Hmgn3, Ilf3, Kcnh7,\n"
    "Malat1, Mbnl2, Pnisr, Prpf39, Rbm34, Rsrc2, Slc25a40, Tia1, Zranb2.\n"
    "Metabolism: Slc5a8, Emb, Bsg, Slc4a4, Crat, Crot, Slc10a2, Acat1, Acat2, Gapdh,\n"
    "Slc13a5, Pdha1, Pdhb, Gpx4, Gpx1, Glul, Adh1, Abcc3, Abhd14b, Pfkm, Cfhr1,\n"
    "Dgkb, Gng12, Elovl5."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 250),
                    genes_text, fontsize=8.5, fontname='Helvetica')

y = 730
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Caveats", fontsize=10, fontname='Helvetica-Bold')
y += 16
caveats = (
    "Preliminary. n=2-3 mice per group. Bars show per-mouse mean ± SEM; "
    "dots = individual mice. Bars/SEM with single replicate (n=1 group from one slide) "
    "have no error bar."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    caveats, fontsize=8, fontname='Helvetica-Oblique')

COVER = BASE / '_subtype_barplot_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()

# --- Merge ---
out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))

# Section page: Excitatory
section = fitz.open()
sp = section.new_page(width=612, height=792)
sp.insert_textbox(fitz.Rect(54, 300, 558, 380),
                  "PART 1 — Excitatory neurons\n"
                  "(7 pages; 122 genes; 20 genes per page)",
                  fontsize=18, fontname='Helvetica-Bold')
SECTION_E = BASE / '_section_excit_tmp.pdf'
section.save(str(SECTION_E))
section.close()
out.insert_pdf(fitz.open(str(SECTION_E)))
out.insert_pdf(fitz.open(str(EXCIT)))

# Section page: Inhibitory
section = fitz.open()
sp = section.new_page(width=612, height=792)
sp.insert_textbox(fitz.Rect(54, 300, 558, 380),
                  "PART 2 — Inhibitory neurons\n"
                  "(7 pages; 122 genes; 20 genes per page)",
                  fontsize=18, fontname='Helvetica-Bold')
SECTION_I = BASE / '_section_inhib_tmp.pdf'
section.save(str(SECTION_I))
section.close()
out.insert_pdf(fitz.open(str(SECTION_I)))
out.insert_pdf(fitz.open(str(INHIB)))

out.save(str(OUT))
out.close()
print(f'Wrote {OUT}')

for tmp in (COVER, SECTION_E, SECTION_I):
    if tmp.exists(): os.remove(tmp)
