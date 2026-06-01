#!/usr/bin/env python3
"""Build a combined PDF for the 3-subtype side-by-side comparison."""
from __future__ import annotations

import os
from datetime import date
from pathlib import Path

import fitz

BASE = Path('/Users/naomi/Desktop/10x/10 Xenium may 2026 data/analysis/current')
DATE = date.today().isoformat()
OUT = BASE / f'NeuronAstrocyte_3SubtypeCompare_Summary_{DATE}.pdf'
SRC = BASE / f'SlidesAB_3subtype_compare_BarPlots_{DATE}.pdf'

cover = fitz.open()
page = cover.new_page(width=612, height=792)

page.insert_textbox(fitz.Rect(54, 60, 558, 110),
                    "Hippocampal Excitatory | Inhibitory | Astrocyte —\n"
                    "side-by-side per-mouse bar plots across all 7 condition groups",
                    fontsize=14, fontname='Helvetica-Bold')
page.insert_textbox(fitz.Rect(54, 112, 558, 128),
                    f"Naomi Kassahun · {DATE} · MEWS Lab · Xenium May 2026 (Slides A + B)",
                    fontsize=10, fontname='Helvetica-Oblique')

y = 150
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Layout", fontsize=12, fontname='Helvetica-Bold')
y += 22
layout = (
    "Each row shows one gene with three small bar-plot panels side-by-side:\n"
    "  Excitatory neurons  |  Inhibitory neurons  |  Astrocytes\n"
    "Each panel has 7 bars (one per condition group). Bars: per-mouse mean.\n"
    "Error bars: SEM across mice. Black dots: individual mouse means.\n"
    "Within each row, the three panels share the same y-axis range so cell-type\n"
    "comparisons are direct.\n"
    "6 genes per page x 40 pages = 236 genes plotted."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 110),
                    layout, fontsize=9, fontname='Helvetica')

y = 268
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Cell-type definitions",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
ct_text = (
    "Excitatory: cells with Neuron lineage score >0 AND Slc17a7/Camk2a/Camk2b/\n"
    "    Satb2/Tbr1/Neurod6 module score >0. ~75,000 cells.\n"
    "Inhibitory: cells with Neuron lineage score >0 AND Gad1/Gad2/Slc32a1/Pvalb/\n"
    "    Sst/Vip/Reln/Lhx6 module score >0. ~32,000 cells.\n"
    "Astrocyte: cells whose top-scoring lineage is the Astrocyte module\n"
    "    (Gfap/Aqp4/Slc1a3/Aldh1l1/S100b/Aldoc/Gja1). ~44,000 cells."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 100),
                    ct_text, fontsize=9, fontname='Helvetica')

y = 370
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Gene list (236 genes plotted; all 95 user-listed genes confirmed present)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
genes_text = (
    "1. Cell-type identity (sanity): Rbfox3 (NeuN), Snap25, Syn1, Map2;\n"
    "   Slc17a7 (VGLUT1), Slc17a6 (VGLUT2), Camk2a, Camk2b, Satb2, Tbr1, Neurod6;\n"
    "   Gad1, Gad2, Slc32a1 (VGAT), Pvalb, Sst, Vip, Reln;\n"
    "   Gfap, Aqp4, Slc1a3 (GLAST), Aldh1l1, S100b, Aldoc, Gja1, Slc1a1 (EAAC1), Slc6a11 (GAT-3).\n"
    "2. MCT family: Slc16a1 (MCT1), Slc16a7 (MCT2), Slc16a3 (MCT4).\n"
    "3. Glucose transporters: Slc2a1 (GLUT1), Slc2a3 (GLUT3), Slc2a4 (GLUT4),\n"
    "   Slc2a2 (GLUT2), Slc2a5 (GLUT5), Slc5a2 (SGLT2).\n"
    "4. Reporter: GFP.\n"
    "5. Alcohol-metabolizing enzymes: Adh1, Adh5, Aldh2, Aldh1a2, Cyp2e1, Cat (catalase).\n"
    "6. Acetate / Ac-CoA / glycolysis-TCA / lactate: Acss1, Acss2, Acly, Eno1, Eno2, Ldha, Ldhb.\n"
    "7. Oxidative stress / detox: Sod1, Sod2, Gpx1, Gpx4, Nfe2l2, Nfe2l1, Hmox1, Keap1,\n"
    "   Sirt3, Sirt5, Gstp1.\n"
    "8. NEW Hormone receptors: Nr3c1 (GR), Nr3c2 (MR), Esr1, Esr2, Ar, Pgr,\n"
    "   Thra (TRalpha), Thrb (TRbeta), Insr (INSR), Igf1r, Igf2r, Lepr, Ghr,\n"
    "   Crhr1, Crhr2, Oxtr, Avpr1b, Mc4r.\n"
    "9. NEW Metabolic / ER stress / UPR / ISR: Atf4 (ATF4), Atf6 (ATF6), Ddit3 (CHOP),\n"
    "   Hspa5 (BiP), Eif2ak3 (PERK), Eif2ak2 (PKR), Eif2ak4 (GCN2), Eif2s1 (eIF2alpha),\n"
    "   Sesn1, Sesn2, Tsc1, Tsc2, Mtor (mTOR), Rheb, Foxo1, Foxo3, Foxo4,\n"
    "   Pmaip1 (NOXA), Bbc3 (PUMA), Trib3.\n"
    "10. NEW Comprehensive IEGs (immediate-early genes):\n"
    "    Fos family: Fos, Fosb, Fosl1, Fosl2.\n"
    "    Jun family: Jun, Junb.\n"
    "    Egr family: Egr1, Egr2, Egr3, Egr4.\n"
    "    Nr4a family: Nr4a1 (Nur77), Nr4a2 (Nurr1), Nr4a3 (NOR1).\n"
    "    Activity-regulated: Arc, Bdnf, Homer1, Npas4, Per1, Per2,\n"
    "      Dusp1, Dusp6, Plk2, Plk3, Nptx1, Nptx2, Trib2, Trib3.\n"
    "    Other IEGs: Atf3, Btg2, Gadd45b, Bhlhe40, Klf10, Ccn1, Rasgrf1.\n"
    "11. NEW Neurotrophic factors + receptors:\n"
    "    Factors: Bdnf, Ntf3 (NT-3), Ntf5 (NT-4/5), Ngf, Gdnf, Cntf, Lif,\n"
    "      Mst1, Mstn, Igf1, Igf2.\n"
    "    Receptors: Ntrk1 (TrkA), Ntrk2 (TrkB), Ntrk3 (TrkC), Ngfr (p75-NTR),\n"
    "      Gfra1, Gfra2, Gfra3, Ret, Lifr.\n"
    "12. User-supplied 95-gene custom list (chromatin/methylation/one-carbon,\n"
    "    synaptic, RNA-binding, metabolism, signaling) - all present in panel."
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 260),
                    genes_text, fontsize=8.5, fontname='Helvetica')

y = 660
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Sample sizes (mice per group)",
                    fontsize=12, fontname='Helvetica-Bold')
y += 22
samples = (
    "  H2O_veh: 2     H2O_MCT1i: 2     EtOH_veh: 2     EtOH_MCT1i: 3\n"
    "  ChronicEtOH: 3     MAT2A_CM: 3     MAT2A_OE: 3"
)
page.insert_textbox(fitz.Rect(54, y, 558, y + 50),
                    samples, fontsize=9, fontname='Courier')

y = 725
page.insert_textbox(fitz.Rect(54, y, 558, y + 18),
                    "Caveats", fontsize=10, fontname='Helvetica-Bold')
y += 16
caveats = (
    "Preliminary. n=2-3 mice per group; bars show per-mouse mean ± SEM with dots = "
    "individual mice. Genes with mostly-zero expression in a given cell type appear as "
    "very small or absent bars — check the per-mouse dots for the actual signal."
)
page.insert_textbox(fitz.Rect(54, y, 558, 770),
                    caveats, fontsize=8, fontname='Helvetica-Oblique')

COVER = BASE / '_compare_cover_tmp.pdf'
cover.save(str(COVER))
cover.close()

out = fitz.open()
out.insert_pdf(fitz.open(str(COVER)))
out.insert_pdf(fitz.open(str(SRC)))
out.save(str(OUT))
out.close()
print(f'Wrote {OUT}')
os.remove(COVER)
