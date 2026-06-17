"""
Translocation v3 — Publication-ready heatmaps (FC vs Naive)
=============================================================
Heatmap A  — 346 proteins, 3 conditions
Heatmap B  — 60 true translocation proteins, 3 conditions
Heatmap C  — 346 proteins, AW only
Heatmap D  — 60 true translocation proteins, AW only
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
STATS_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_v3_stats.xlsx')

# ── Colormap ──────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'],
    N=512,
)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

# ── Fixed layout constants ────────────────────────────────────────────────────
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

# ── Core heatmap function ─────────────────────────────────────────────────────
def make_heatmap(df, title, out_pdf, out_png,
                 cond_labels, panels_fn,
                 inches_per_prot=0.0072,
                 row_h=0.50):
    """
    inches_per_prot : width per protein column (increase for small n to avoid blocky look)
    row_h           : height per condition row in inches
    """
    df     = df.sort_values('Chrom_FC_AW', ascending=False).reset_index(drop=True)
    n_prot = len(df)
    n_cond = len(cond_labels)
    band_h = n_cond * row_h

    hm_w   = n_prot * inches_per_prot
    panels = panels_fn(df)

    FIG_W = LEFT_IN + hm_w + RIGHT_IN
    FIG_H = (TOP_PAD
             + len(panels) * band_h
             + (len(panels) - 1) * SECTION_GAP
             + BOT_PAD)

    fig      = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    im_ref   = None
    y_cursor = FIG_H - TOP_PAD

    for idx, (label, mat) in enumerate(panels):
        y_band_bot = y_cursor - band_h

        ax = fig.add_axes([LEFT_IN / FIG_W,
                           y_band_bot / FIG_H,
                           hm_w / FIG_W,
                           band_h / FIG_H],
                          facecolor=COMP_BG[label])
        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                       interpolation='nearest')
        if im_ref is None:
            im_ref = im

        # Row dividers
        for r in range(1, n_cond):
            ax.axhline(r - 0.5, color='white', linewidth=1.2, zorder=3)

        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        # Condition labels on right
        for r, cname in enumerate(cond_labels):
            y_frac = (y_band_bot + band_h - (r + 0.5) * row_h) / FIG_H
            fig.text((LEFT_IN + hm_w + 0.12) / FIG_W, y_frac,
                     cname, ha='left', va='center',
                     fontsize=8, fontfamily='Arial', color='#222222')

        # Compartment label on left
        mid_y = (y_band_bot + band_h / 2) / FIG_H
        fig.text((LEFT_IN - 0.12) / FIG_W, mid_y,
                 label, ha='right', va='center',
                 fontsize=10, fontweight='bold',
                 fontfamily='Arial', color='#111111')

        # n= below label
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

    # Footnote
    fig.text(0.5, 0.01,
             'AW-M-3 excluded.  Sorted by Chromatin AW FC (high to low).  '
             'Interaction test: delta_CH vs delta_SN, Welch t-test BH FDR p_adj < 0.05.',
             ha='center', va='bottom', fontsize=5.5,
             fontfamily='Arial', color='#888888')

    # Title
    fig.text(0.5, (FIG_H - TOP_PAD * 0.35) / FIG_H,
             title, ha='center', va='top',
             fontsize=9, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_png}')


# ── Panel builders ────────────────────────────────────────────────────────────
def panels_3cond(df):
    return [
        ('Chromatin',       np.vstack([df['Chrom_FC_Intox'].values,
                                       df['Chrom_FC_AW'].values,
                                       df['Chrom_FC_PA'].values])),
        ('Soluble Nuclear', np.vstack([df['SN_FC_Intox'].values,
                                       df['SN_FC_AW'].values,
                                       df['SN_FC_PA'].values])),
    ]

def panels_AW_only(df):
    return [
        ('Chromatin',       df['Chrom_FC_AW'].values[np.newaxis, :]),
        ('Soluble Nuclear', df['SN_FC_AW'].values[np.newaxis, :]),
    ]


# ── Load data ─────────────────────────────────────────────────────────────────
all346 = pd.read_excel(STATS_XLSX, sheet_name='Into_Chromatin')
true60 = pd.read_excel(STATS_XLSX, sheet_name='G1_CH_up_SN_down_n60')

LABELS_3 = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
LABELS_AW = ['Acute Withdrawal']

# ── A: 346 proteins, 3 conditions ────────────────────────────────────────────
make_heatmap(
    all346,
    title           = 'Proteins moved into Chromatin during Acute Withdrawal  (n = 346)',
    out_pdf         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_346_pub.pdf'),
    out_png         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_346_pub.png'),
    cond_labels     = LABELS_3,
    panels_fn       = panels_3cond,
    inches_per_prot = 0.0072,
    row_h           = 0.50,
)

# ── B: 60 true translocation, 3 conditions ───────────────────────────────────
make_heatmap(
    true60,
    title           = 'True translocation into Chromatin during Acute Withdrawal  (n = 60)',
    out_pdf         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_60_pub.pdf'),
    out_png         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_60_pub.png'),
    cond_labels     = LABELS_3,
    panels_fn       = panels_3cond,
    inches_per_prot = 0.045,   # wider columns so 60 proteins fill naturally
    row_h           = 0.50,
)

# ── C: 346 proteins, AW only ──────────────────────────────────────────────────
make_heatmap(
    all346,
    title           = 'Proteins moved into Chromatin — Acute Withdrawal only  (n = 346)',
    out_pdf         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_346_AW_pub.pdf'),
    out_png         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_346_AW_pub.png'),
    cond_labels     = LABELS_AW,
    panels_fn       = panels_AW_only,
    inches_per_prot = 0.0072,
    row_h           = 0.30,    # shorter rows since only 1 condition
)

# ── D: 60 true translocation, AW only ────────────────────────────────────────
make_heatmap(
    true60,
    title           = 'True translocation into Chromatin — Acute Withdrawal only  (n = 60)',
    out_pdf         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_60_AW_pub.pdf'),
    out_png         = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_60_AW_pub.png'),
    cond_labels     = LABELS_AW,
    panels_fn       = panels_AW_only,
    inches_per_prot = 0.045,
    row_h           = 0.30,
)
