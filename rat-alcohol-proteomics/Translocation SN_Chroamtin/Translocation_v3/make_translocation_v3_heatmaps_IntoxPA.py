"""
Translocation v3 — Heatmaps for Intoxication and PA true translocation
=======================================================================
Same aesthetic as make_translocation_v3_heatmaps_pub.py (AW version).

Produces (all 3 conditions shown so temporal context is visible):
  Translocation_Intox_heatmap_all_pub.pdf/.png      — all 76 into Chromatin
  Translocation_Intox_heatmap_true_pub.pdf/.png     — 13 true translocation
  Translocation_PA_heatmap_all_pub.pdf/.png         — all 7 into Chromatin
  Translocation_PA_heatmap_true_pub.pdf/.png        — 3 true translocation
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STATS_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_IntoxPA_stats.xlsx')

# ── Colormap ──────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0','#2B7FD4','#6AAEE0','#AECFE8','#DDEEF8',
     '#FFFFFF',
     '#FDD8E7','#F5A0BC','#EE5F8B','#E8305A','#C01E42'],
    N=512,
)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

# ── Layout ────────────────────────────────────────────────────────────────────
ROW_H       = 0.50
N_COND      = 3
BAND_H      = N_COND * ROW_H
LEFT_IN     = 2.40
RIGHT_IN    = 2.70
TOP_PAD     = 0.55
BOT_PAD     = 1.20
SECTION_GAP = 0.32
CBAR_W_IN   = 2.8
CBAR_H_IN   = 0.14
CBAR_BOT    = 0.48

COMP_BG = {
    'Chromatin':       '#FBF9F0',
    'Soluble Nuclear': '#F0F6FB',
}
COND_LABELS = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

# ── Core draw function ────────────────────────────────────────────────────────
def make_heatmap(df, title, sort_col, out_pdf, out_png, inches_per_prot=0.0072):
    """
    sort_col : which Chrom_FC column to sort by (Chrom_FC_Intox or Chrom_FC_PA)
    """
    df     = df.sort_values(sort_col, ascending=False).reset_index(drop=True)
    n_prot = len(df)

    panels = [
        ('Chromatin',       np.vstack([df['Chrom_FC_Intox'].values,
                                       df['Chrom_FC_AW'].values,
                                       df['Chrom_FC_PA'].values])),
        ('Soluble Nuclear', np.vstack([df['SN_FC_Intox'].values,
                                       df['SN_FC_AW'].values,
                                       df['SN_FC_PA'].values])),
    ]

    hm_w  = max(n_prot * inches_per_prot, 2.0)   # minimum 2 inches for small n
    FIG_W = LEFT_IN + hm_w + RIGHT_IN
    FIG_H = (TOP_PAD
             + len(panels) * BAND_H
             + (len(panels) - 1) * SECTION_GAP
             + BOT_PAD)

    fig      = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    im_ref   = None
    y_cursor = FIG_H - TOP_PAD

    for label, mat in panels:
        y_band_bot = y_cursor - BAND_H

        ax = fig.add_axes([LEFT_IN / FIG_W, y_band_bot / FIG_H,
                           hm_w / FIG_W,    BAND_H / FIG_H],
                          facecolor=COMP_BG[label])

        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                       interpolation='nearest')
        if im_ref is None:
            im_ref = im

        for r in range(1, N_COND):
            ax.axhline(r - 0.5, color='white', linewidth=1.2, zorder=3)

        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        # Condition labels
        for r, cname in enumerate(COND_LABELS):
            y_frac = (y_band_bot + BAND_H - (r + 0.5) * ROW_H) / FIG_H
            fig.text((LEFT_IN + hm_w + 0.12) / FIG_W, y_frac,
                     cname, ha='left', va='center',
                     fontsize=8, fontfamily='Arial', color='#222222')

        # Compartment label
        mid_y = (y_band_bot + BAND_H / 2) / FIG_H
        fig.text((LEFT_IN - 0.12) / FIG_W, mid_y,
                 label, ha='right', va='center',
                 fontsize=10, fontweight='bold',
                 fontfamily='Arial', color='#111111')

        fig.text((LEFT_IN - 0.12) / FIG_W, mid_y - 0.055,
                 f'n = {n_prot}', ha='right', va='top',
                 fontsize=7.5, fontfamily='Arial',
                 color='#555555', style='italic')

        y_cursor = y_band_bot - SECTION_GAP

    # Colorbar
    cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
    ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H,
                           CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
    cb = fig.colorbar(im_ref, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
    cb.ax.tick_params(labelsize=7)
    ax_cb.set_xlabel('Log2 Fold Change vs Naive',
                     fontsize=8, fontfamily='Arial', color='#444444', labelpad=4)

    fig.text(0.5, 0.01,
             'AW-M-3 excluded.  Sorted by Chromatin FC of condition (high to low).  '
             'Interaction test: delta_CH vs delta_SN, Welch t-test BH FDR p_adj < 0.05.',
             ha='center', va='bottom', fontsize=5.5,
             fontfamily='Arial', color='#888888')

    fig.text(0.5, (FIG_H - TOP_PAD * 0.35) / FIG_H,
             title, ha='center', va='top',
             fontsize=9, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_png}')


# ── Load data ─────────────────────────────────────────────────────────────────
intox_all  = pd.read_excel(STATS_XLSX, sheet_name='Into_Chromatin_Intox')
intox_true = pd.read_excel(STATS_XLSX, sheet_name='G1_true_translocation_Intox')
pa_all     = pd.read_excel(STATS_XLSX, sheet_name='Into_Chromatin_PA')
pa_true    = pd.read_excel(STATS_XLSX, sheet_name='G1_true_translocation_PA')

print(f'Intox into Chromatin: {len(intox_all)}  |  True translocation: {len(intox_true)}')
print(f'PA into Chromatin:    {len(pa_all)}      |  True translocation: {len(pa_true)}')

# ── Intox: all 76 ─────────────────────────────────────────────────────────────
make_heatmap(
    intox_all,
    title          = f'Proteins moved into Chromatin during Intoxication  (n = {len(intox_all)})',
    sort_col       = 'Chrom_FC_Intox',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_all_pub.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_all_pub.png'),
    inches_per_prot= 0.0072,
)

# ── Intox: true translocation (13) ────────────────────────────────────────────
make_heatmap(
    intox_true,
    title          = f'True translocation into Chromatin during Intoxication  (n = {len(intox_true)})',
    sort_col       = 'Chrom_FC_Intox',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_true_pub.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_true_pub.png'),
    inches_per_prot= 0.10,   # wider columns for small n
)

# ── PA: all 7 ─────────────────────────────────────────────────────────────────
make_heatmap(
    pa_all,
    title          = f'Proteins moved into Chromatin during Protracted Abstinence  (n = {len(pa_all)})',
    sort_col       = 'Chrom_FC_PA',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_all_pub.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_all_pub.png'),
    inches_per_prot= 0.10,
)

# ── PA: true translocation (3) ────────────────────────────────────────────────
make_heatmap(
    pa_true,
    title          = f'True translocation into Chromatin during Protracted Abstinence  (n = {len(pa_true)})',
    sort_col       = 'Chrom_FC_PA',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_true_pub.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_true_pub.png'),
    inches_per_prot= 0.30,   # very wide for only 3 proteins
)
