"""
Translocation v3 — Into Soluble Nuclear (reverse direction)
============================================================
Reads the already-computed Into_SN sheets from existing stats Excel files,
applies the true translocation filter (SN up + Chrom_FC < -0.5 for the
relevant condition), generates heatmaps for all 3 conditions, and saves a
summary Excel file and a 3×2 summary figure.

True translocation out of Chromatin into SN:
  - Significant interaction score favouring SN (already in Into_SN sheets)
  - Chrom_FC_<cond> < -0.5  (chromatin depleted = protein moved OUT of chromatin)

Outputs (all in same folder as this script):
  Translocation_SN_stats.xlsx                     — 6 sheets, all conditions
  Translocation_SN_Intox_heatmap_all.pdf/.png      — 41 proteins
  Translocation_SN_Intox_heatmap_true.pdf/.png     — 34 proteins
  Translocation_SN_AW_heatmap_all.pdf/.png         — 362 proteins
  Translocation_SN_AW_heatmap_true.pdf/.png        — 355 proteins
  Translocation_SN_PA_heatmap_all.pdf/.png         — 7 proteins
  Translocation_SN_PA_heatmap_true.pdf/.png        — 7 proteins
  Translocation_SN_Summary_AllConditions.pdf/.png  — 3×2 summary figure
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
import matplotlib.image as mpimg

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
STATS_AW    = os.path.join(SCRIPT_DIR, 'Translocation_v3_stats.xlsx')
STATS_IP    = os.path.join(SCRIPT_DIR, 'Translocation_IntoxPA_stats.xlsx')
OUT_XLSX    = os.path.join(SCRIPT_DIR, 'Translocation_SN_stats.xlsx')

# ── Colormap (identical to v3) ────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0','#2B7FD4','#6AAEE0','#AECFE8','#DDEEF8',
     '#FFFFFF',
     '#FDD8E7','#F5A0BC','#EE5F8B','#E8305A','#C01E42'],
    N=512,
)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

# ── Heatmap layout (same as IntoxPA heatmaps) ─────────────────────────────────
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

# ── Load and build subsets ────────────────────────────────────────────────────
print('Loading Into_SN sheets...')

# AW
sn_aw_all  = pd.read_excel(STATS_AW, sheet_name='Into_SN')
# True translocation = SN actually went UP + Chromatin actually went DOWN
# Mirror of Into-Chromatin true: CH_FC > 0 AND SN_FC < -0.5
# Into-SN true:                  SN_FC > 0 AND Chrom_FC < -0.5
sn_aw_true = sn_aw_all[
    (sn_aw_all['SN_FC_AW'] > 0) & (sn_aw_all['Chrom_FC_AW'] < -0.5)
].copy()

# Intox
sn_i_all   = pd.read_excel(STATS_IP, sheet_name='Into_SN_Intox')
sn_i_true  = sn_i_all[
    (sn_i_all['SN_FC_Intox'] > 0) & (sn_i_all['Chrom_FC_Intox'] < -0.5)
].copy()

# PA
sn_p_all   = pd.read_excel(STATS_IP, sheet_name='Into_SN_PA')
sn_p_true  = sn_p_all[
    (sn_p_all['SN_FC_PA'] > 0) & (sn_p_all['Chrom_FC_PA'] < -0.5)
].copy()

print(f'Intox  — all: {len(sn_i_all):>4}  |  true translocation (SN↑ & CH↓): {len(sn_i_true)}')
print(f'AW     — all: {len(sn_aw_all):>4}  |  true translocation (SN↑ & CH↓): {len(sn_aw_true)}')
print(f'PA     — all: {len(sn_p_all):>4}  |  true translocation (SN↑ & CH↓): {len(sn_p_true)}')

# ── Save Excel ────────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    sn_i_all.to_excel(writer,   sheet_name='Into_SN_Intox',              index=False)
    sn_i_true.to_excel(writer,  sheet_name='G1_true_translocation_Intox',index=False)
    sn_aw_all.to_excel(writer,  sheet_name='Into_SN_AW',                 index=False)
    sn_aw_true.to_excel(writer, sheet_name='G1_true_translocation_AW',   index=False)
    sn_p_all.to_excel(writer,   sheet_name='Into_SN_PA',                 index=False)
    sn_p_true.to_excel(writer,  sheet_name='G1_true_translocation_PA',   index=False)
print(f'Saved: {OUT_XLSX}')

# ── Core heatmap function ─────────────────────────────────────────────────────
def make_heatmap(df, title, sort_col, out_pdf, out_png, inches_per_prot=0.0072):
    """
    sort_col : SN_FC column for the primary condition (sort descending = most into SN at left)
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

    hm_w  = max(n_prot * inches_per_prot, 2.0)
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
             'AW-M-3 excluded.  Sorted by SN FC of condition (high to low).  '
             'Interaction test: delta_SN vs delta_CH, Welch t-test BH FDR p_adj < 0.05.',
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


# ── Generate all 6 heatmaps ───────────────────────────────────────────────────

# Intox — all (41)
make_heatmap(
    sn_i_all,
    title          = f'Proteins moved into Soluble Nuclear during Intoxication  (n = {len(sn_i_all)})',
    sort_col       = 'SN_FC_Intox',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_all.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_all.png'),
    inches_per_prot= 0.10,
)

# Intox — true (13)
make_heatmap(
    sn_i_true,
    title          = f'True translocation into Soluble Nuclear during Intoxication  (n = {len(sn_i_true)})',
    sort_col       = 'SN_FC_Intox',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_true.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_true.png'),
    inches_per_prot= 0.20,
)

# AW — all (362)
make_heatmap(
    sn_aw_all,
    title          = f'Proteins moved into Soluble Nuclear during Acute Withdrawal  (n = {len(sn_aw_all)})',
    sort_col       = 'SN_FC_AW',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_all.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_all.png'),
    inches_per_prot= 0.0072,
)

# AW — true (131)
make_heatmap(
    sn_aw_true,
    title          = f'True translocation into Soluble Nuclear during Acute Withdrawal  (n = {len(sn_aw_true)})',
    sort_col       = 'SN_FC_AW',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_true.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_true.png'),
    inches_per_prot= 0.035,
)

# PA — all (7)
make_heatmap(
    sn_p_all,
    title          = f'Proteins moved into Soluble Nuclear during Protracted Abstinence  (n = {len(sn_p_all)})',
    sort_col       = 'SN_FC_PA',
    out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_all.pdf'),
    out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_all.png'),
    inches_per_prot= 0.20,
)

# PA — true (0 proteins — no genuine Chrom→SN translocation during abstinence)
# Generate a blank placeholder PNG instead of a heatmap
if len(sn_p_true) == 0:
    fig_blank = plt.figure(figsize=(6, 4), facecolor='white')
    fig_blank.text(0.5, 0.55,
                   'No true translocation detected',
                   ha='center', va='center', fontsize=11,
                   fontfamily='Arial', color='#888888', style='italic')
    fig_blank.text(0.5, 0.42,
                   'n = 0  (SN↑ & Chrom↓ during Protracted Abstinence)',
                   ha='center', va='center', fontsize=8,
                   fontfamily='Arial', color='#AAAAAA')
    _pa_true_png = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_true.png')
    _pa_true_pdf = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_true.pdf')
    with PdfPages(_pa_true_pdf) as pdf:
        pdf.savefig(fig_blank, dpi=200, bbox_inches='tight', facecolor='white')
    fig_blank.savefig(_pa_true_png, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig_blank)
    print(f'Saved placeholder: {_pa_true_png}')
else:
    make_heatmap(
        sn_p_true,
        title          = f'True translocation into Soluble Nuclear during Protracted Abstinence  (n = {len(sn_p_true)})',
        sort_col       = 'SN_FC_PA',
        out_pdf        = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_true.pdf'),
        out_png        = os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_true.png'),
        inches_per_prot= 0.30,
    )

# ── 3×2 Summary figure ────────────────────────────────────────────────────────
print('\nBuilding summary figure...')

IMAGES = [
    [os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_all.png'),
     os.path.join(SCRIPT_DIR, 'Translocation_SN_Intox_heatmap_true.png')],
    [os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_all.png'),
     os.path.join(SCRIPT_DIR, 'Translocation_SN_AW_heatmap_true.png')],
    [os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_all.png'),
     os.path.join(SCRIPT_DIR, 'Translocation_SN_PA_heatmap_true.png')],
]

ROW_LABELS  = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
COL_LABELS  = ['All proteins into Soluble Nuclear',
               'True Translocation\n(SN ↑,  Chromatin ↓)']
ROW_COLORS  = ['#FFF3E0', '#E8F4FB', '#F0FFF0']

FIG_W  = 18.0
FIG_H  = 22.0
LEFT_LABEL_W = 1.0
TOP_LABEL_H  = 0.55
GAP_X        = 0.30
GAP_Y        = 0.45

cell_w = (FIG_W - LEFT_LABEL_W - GAP_X) / 2
cell_h = (FIG_H - TOP_LABEL_H  - 2 * GAP_Y) / 3

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

# Column headers
for col, col_label in enumerate(COL_LABELS):
    x_center = (LEFT_LABEL_W + col * (cell_w + GAP_X) + cell_w / 2) / FIG_W
    y_pos    = (FIG_H - TOP_LABEL_H * 0.45) / FIG_H
    fig.text(x_center, y_pos, col_label,
             ha='center', va='center',
             fontsize=13, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')

for row in range(3):
    y_top = FIG_H - TOP_LABEL_H - row * (cell_h + GAP_Y)
    y_bot = y_top - cell_h

    # Row background
    strip = FancyBboxPatch(
        (0, y_bot / FIG_H), 1.0, cell_h / FIG_H,
        boxstyle='square,pad=0', facecolor=ROW_COLORS[row],
        edgecolor='none', transform=fig.transFigure, zorder=0, clip_on=False,
    )
    fig.add_artist(strip)

    # Row label
    fig.text(
        (LEFT_LABEL_W * 0.42) / FIG_W,
        (y_bot + cell_h / 2) / FIG_H,
        ROW_LABELS[row],
        ha='center', va='center',
        fontsize=12, fontweight='bold',
        fontfamily='Arial', color='#1A1A1A', rotation=90,
    )

    # Vertical separator
    fig.add_artist(plt.matplotlib.lines.Line2D(
        [LEFT_LABEL_W / FIG_W, LEFT_LABEL_W / FIG_W],
        [y_bot / FIG_H, y_top / FIG_H],
        transform=fig.transFigure,
        color='#CCCCCC', linewidth=0.8, zorder=1,
    ))

    for col in range(2):
        x_left = LEFT_LABEL_W + col * (cell_w + GAP_X)

        # Border
        border = FancyBboxPatch(
            (x_left / FIG_W, y_bot / FIG_H),
            cell_w / FIG_W, cell_h / FIG_H,
            boxstyle='round,pad=0.003',
            facecolor='none', edgecolor='#DDDDDD', linewidth=0.8,
            transform=fig.transFigure, zorder=3, clip_on=False,
        )
        fig.add_artist(border)

        # Image
        img_path = IMAGES[row][col]
        if os.path.exists(img_path):
            img = mpimg.imread(img_path)
            ax  = fig.add_axes([x_left / FIG_W, y_bot / FIG_H,
                                 cell_w / FIG_W, cell_h / FIG_H])
            ax.imshow(img, aspect='auto')
            ax.axis('off')
        else:
            ax = fig.add_axes([x_left / FIG_W, y_bot / FIG_H,
                                cell_w / FIG_W, cell_h / FIG_H])
            ax.set_facecolor('#F5F5F5')
            ax.text(0.5, 0.5, 'File not found', ha='center', va='center',
                    fontsize=10, color='#999999')
            ax.axis('off')

# Horizontal dividers
for row in range(1, 3):
    y = (FIG_H - TOP_LABEL_H - row * (cell_h + GAP_Y) + GAP_Y / 2) / FIG_H
    fig.add_artist(plt.matplotlib.lines.Line2D(
        [0, 1], [y, y],
        transform=fig.transFigure,
        color='#CCCCCC', linewidth=0.8, zorder=2,
    ))

# Title
fig.text(0.5, (FIG_H - 0.15) / FIG_H,
         'Protein Translocation into Soluble Nuclear — All Conditions',
         ha='center', va='top',
         fontsize=15, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')

# Footnote
fig.text(0.5, 0.005,
         'Delta-delta interaction test (delta_SN vs delta_CH, Welch t-test, BH FDR p_adj < 0.05).  '
         'AW-M-3 excluded.  Chromatin: Keep + Review only.  Color = Log2 FC vs Naive.',
         ha='center', va='bottom',
         fontsize=7, fontfamily='Arial', color='#888888')

out_sum_pdf = os.path.join(SCRIPT_DIR, 'Translocation_SN_Summary_AllConditions.pdf')
out_sum_png = os.path.join(SCRIPT_DIR, 'Translocation_SN_Summary_AllConditions.png')

with PdfPages(out_sum_pdf) as pdf:
    pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(out_sum_png, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_sum_png}')
print('\nAll done.')
