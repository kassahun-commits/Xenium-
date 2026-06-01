#!/usr/bin/env python3
"""Cover + 51-page grouped-celltype bar plot grid."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()
OUT = BASE / f'NeuronAstrocyte_GroupedCellType_BarPlot_Summary_{DATE}.pdf'
SRC = BASE / f'SlidesAB_grouped_celltype_BarPlots_{DATE}.pdf'

cover = fitz.open()
page = cover.new_page(width=612, height=792)
page.insert_textbox(fitz.Rect(54, 60, 558, 115),
                    "Hippocampal Excitatory | Inhibitory | Astrocyte —\n"
                    "grouped-by-condition bar plots\n"
                    "(3 cell types adjacent within each condition group)",
                    fontsize=13, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 116, 558, 132),
                    f"Naomi Kassahun · {DATE} · MEWS Lab · Xenium May 2026 (Slides A + B)",
                    fontsize=10, fontname='Helvetica-Oblique')

y = 155
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Layout", fontsize=12, fontname='Helvetica-Bold')
y += 22
layout = (
    "Each panel = one gene with all 21 bars.\n"
    "Along x: 7 condition-group clusters (H2O_veh, H2O_MCT1i, EtOH_veh,\n"
    "EtOH_MCT1i, ChronicEtOH, MAT2A_CM, MAT2A_OE).\n"
    "Within each cluster: 3 bars side-by-side colored by cell type —\n"
    "  blue = Excitatory, purple = Inhibitory, red = Astrocyte.\n"
    "Bar height = mean of per-mouse means. Error bars = SEM. Dots = individual mice.\n"
    "Light vertical lines separate the 7 condition clusters.\n"
    "\n"
    "6 panels per page (2 columns × 3 rows). 51 plot pages total.\n"
    "Genes are GROUPED BY FUNCTION — one or more pages per category, with\n"
    "the category name shown in the page title."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 195),
                    layout, fontsize=9, fontname='Helvetica')

y = 360
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Functional categories (26 total)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
cats_text = (
    "1.  Neuronal identity (sanity)\n"
    "2.  Astrocyte identity (sanity)\n"
    "3.  MCT family — MCT1, MCT2, MCT4\n"
    "4.  Glucose transporters — GLUT1, GLUT2, GLUT3, GLUT4, GLUT5, SGLT2\n"
    "5.  Reporter — GFP\n"
    "6.  Alcohol-metabolizing enzymes — Adh1, Adh5, Aldh2, Aldh1a2, Cyp2e1, Cat\n"
    "7.  Acetate / Ac-CoA / lactate — Acss1, Acss2, Acly, Ldha, Ldhb\n"
    "8.  Glycolysis / TCA / energy — Gapdh, Pgk1, Pkm, Hk1/2, Pfkm, Aldoa, Eno1/2,\n"
    "     Idh1/2, Cs, Pdha1, Pdhb\n"
    "9.  One-carbon / SAM / methylation metabolism — Mat1a/2a/2b, Ahcy, Mthfd1/2,\n"
    "     Mtr, Mtrr, Pemt, Phgdh, Psph, Shmt2, Tyms, Mtap, Gnmt, Chdh\n"
    "10. Lipid / fatty-acid metabolism — Acat1, Acat2, Crat, Crot, Elovl5\n"
    "11. Solute transporters (misc) — Slc4a4, Slc5a8, Slc10a2, Slc13a5, Slc25a32/40,\n"
    "     Bsg, Emb, Abcc3\n"
    "12. Oxidative stress / detox — Sod1/2, Gpx1/4, Nfe2l1/2, Keap1, Hmox1,\n"
    "     Sirt3/5, Gstp1\n"
    "13. Metabolic / ER stress / UPR / ISR — ATF4/6, CHOP, BiP, PERK, PKR, GCN2,\n"
    "     eIF2α, Sestrin1/2, TSC1/2, mTOR, Rheb, FoxO1/3/4, NOXA, PUMA, Trib3\n"
    "14. Hormone receptors — GR, MR, ERα, ERβ, AR, PR, TRα, TRβ, INSR, IGF1R,\n"
    "     IGF2R, LEPR, GHR, CRHR1/2, OXTR, V1bR, MC4R\n"
    "15. IEGs — Fos family (Fos, Fosb, Fosl1, Fosl2)\n"
    "16. IEGs — Jun family (Jun, Junb)\n"
    "17. IEGs — Egr family (Egr1, Egr2, Egr3, Egr4)\n"
    "18. IEGs — Nr4a family (Nur77, Nurr1, NOR1)\n"
    "19. IEGs — activity-regulated (Arc, Bdnf, Homer1, Npas4, Per1/2, Dusp1/6,\n"
    "     Plk2/3, Nptx1/2, Trib2)\n"
    "20. IEGs — other / activity-related (Atf3, Btg2, Gadd45b, Bhlhe40, Klf10,\n"
    "     Ccn1, Rasgrf1)\n"
    "21. Neurotrophic factors — NT-3, NT-4/5, NGF, GDNF, CNTF, LIF, MST1, MSTN,\n"
    "     IGF1, IGF2\n"
    "22. Neurotrophic receptors — TrkA, TrkB, TrkC, p75NTR, Gfra1-3, Ret, LIFR\n"
    "23. Chromatin / histone modifiers — Brd9, Chd3/8, Dpf1, Ino80, Kat7, Kmt5a,\n"
    "     Mbd1/4, Nsd2/3, Prdm2/8, Prmt2/6/8/9, Setdb2, Ss18, Suv39h2\n"
    "24. Synaptic / structural — Grik5, Rgs14, Pcp4, Cacng5, Cacnb4, Dlg2,\n"
    "     Lrrtm2/4, Pclo, Unc13a, Tagln3\n"
    "25. RNA-binding / nuclear / splicing — Celf2, Ccnt2, Cdk2ap1, Hmgb3, Hmgn3,\n"
    "     Ilf3, Kcnh7, Malat1, Mbnl2, Pnisr, Prpf39, Rbm34, Rsrc2, Tia1, Zranb2\n"
    "26. Other / signaling — Cfhr1, Abhd14b, Dgkb, Gng12"
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 380),
                    cats_text, fontsize=8, fontname='Helvetica')

y = 770
page.insert_textbox(fitz.Rect(54, y, 558, 790),
                    "Sample sizes: H2O_veh n=2, H2O_MCT1i n=2, EtOH_veh n=2, "
                    "EtOH_MCT1i n=3, ChronicEtOH n=3, MAT2A_CM n=3, MAT2A_OE n=3.",
                    fontsize=8, fontname='Helvetica-Oblique')

COVER = BASE / '_grouped_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()

out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))
out.insert_pdf(fitz.open(str(SRC)))
out.save(str(OUT))
out.close()
print(f'Wrote {OUT}')
os.remove(COVER)
