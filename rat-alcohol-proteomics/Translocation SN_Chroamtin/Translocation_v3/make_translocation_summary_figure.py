"""
Translocation v3 — Summary figure combining all heatmaps
=========================================================
Layout: 3 rows (Intox, AW, PA) × 2 columns (All into Chromatin | True Translocation)
Output: Translocation_Summary_AllConditions.pdf/.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
import matplotlib.image as mpimg

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Images to load ────────────────────────────────────────────────────────────
# [row][col]: row = condition, col = 0 (all) or 1 (true translocation)
IMAGES = [
    [  # Intoxication
        os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_all_pub.png'),
        os.path.join(SCRIPT_DIR, 'Translocation_Intox_heatmap_true_pub.png'),
    ],
    [  # Acute Withdrawal
        os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_346_pub.png'),
        os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_60_pub.png'),
    ],
    [  # Protracted Abstinence
        os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_all_pub.png'),
        os.path.join(SCRIPT_DIR, 'Translocation_PA_heatmap_true_pub.png'),
    ],
]

N_LABELS = [('n = 76', 'n = 13'),
            ('n = 346', 'n = 60'),
            ('n = 7',   'n = 3')]

ROW_LABELS = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
COL_LABELS = ['All proteins into Chromatin', 'True Translocation\n(Chromatin ↑,  Soluble Nuclear ↓)']

ROW_COLORS = ['#FFF3E0', '#E8F4FB', '#F0FFF0']   # warm / blue / green tint per condition

# ── Figure layout ─────────────────────────────────────────────────────────────
FIG_W  = 18.0
FIG_H  = 22.0

LEFT_LABEL_W = 1.0    # inches for row label strip
TOP_LABEL_H  = 0.55   # inches for column header
GAP_X        = 0.30   # horizontal gap between columns
GAP_Y        = 0.45   # vertical gap between rows

cell_w = (FIG_W - LEFT_LABEL_W - GAP_X) / 2
cell_h = (FIG_H - TOP_LABEL_H  - 2 * GAP_Y) / 3

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

# ── Column headers ────────────────────────────────────────────────────────────
for col, col_label in enumerate(COL_LABELS):
    x_center = (LEFT_LABEL_W + col * (cell_w + GAP_X) + cell_w / 2) / FIG_W
    y_pos    = (FIG_H - TOP_LABEL_H * 0.45) / FIG_H
    fig.text(x_center, y_pos, col_label,
             ha='center', va='center',
             fontsize=13, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')

# ── Row label strip + images ──────────────────────────────────────────────────
for row in range(3):
    y_top = FIG_H - TOP_LABEL_H - row * (cell_h + GAP_Y)
    y_bot = y_top - cell_h

    # Row background strip (subtle tint across full width)
    strip = FancyBboxPatch(
        (0, y_bot / FIG_H),
        1.0, cell_h / FIG_H,
        boxstyle='square,pad=0',
        facecolor=ROW_COLORS[row],
        edgecolor='none',
        transform=fig.transFigure,
        zorder=0, clip_on=False,
    )
    fig.add_artist(strip)

    # Row label (vertical, centred)
    fig.text(
        (LEFT_LABEL_W * 0.42) / FIG_W,
        (y_bot + cell_h / 2) / FIG_H,
        ROW_LABELS[row],
        ha='center', va='center',
        fontsize=12, fontweight='bold',
        fontfamily='Arial', color='#1A1A1A',
        rotation=90,
    )

    # Thin vertical separator line after row label
    fig.add_artist(plt.matplotlib.lines.Line2D(
        [LEFT_LABEL_W / FIG_W, LEFT_LABEL_W / FIG_W],
        [y_bot / FIG_H, y_top / FIG_H],
        transform=fig.transFigure,
        color='#CCCCCC', linewidth=0.8, zorder=1,
    ))

    for col in range(2):
        x_left = LEFT_LABEL_W + col * (cell_w + GAP_X)

        # Thin border around cell (facecolor='none' so it doesn't cover the image)
        border = FancyBboxPatch(
            (x_left / FIG_W, y_bot / FIG_H),
            cell_w / FIG_W, cell_h / FIG_H,
            boxstyle='round,pad=0.003',
            facecolor='none',
            edgecolor='#DDDDDD', linewidth=0.8,
            transform=fig.transFigure,
            zorder=3, clip_on=False,
        )
        fig.add_artist(border)

        # Load and display image
        img_path = IMAGES[row][col]
        if os.path.exists(img_path):
            img = mpimg.imread(img_path)
            ax = fig.add_axes([
                x_left / FIG_W,
                y_bot  / FIG_H,
                cell_w / FIG_W,
                cell_h / FIG_H,
            ])
            ax.imshow(img, aspect='auto')
            ax.axis('off')
        else:
            # Placeholder if file missing
            ax = fig.add_axes([
                x_left / FIG_W, y_bot / FIG_H,
                cell_w / FIG_W, cell_h / FIG_H,
            ])
            ax.set_facecolor('#F5F5F5')
            ax.text(0.5, 0.5, 'File not found', ha='center', va='center',
                    fontsize=10, color='#999999')
            ax.axis('off')

# ── Horizontal dividers between rows ─────────────────────────────────────────
for row in range(1, 3):
    y = (FIG_H - TOP_LABEL_H - row * (cell_h + GAP_Y) + GAP_Y / 2) / FIG_H
    fig.add_artist(plt.matplotlib.lines.Line2D(
        [0, 1], [y, y],
        transform=fig.transFigure,
        color='#CCCCCC', linewidth=0.8, zorder=2,
    ))

# ── Main title ────────────────────────────────────────────────────────────────
fig.text(0.5, (FIG_H - 0.15) / FIG_H,
         'Protein Translocation into Chromatin — All Conditions',
         ha='center', va='top',
         fontsize=15, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')

# ── Footnote ──────────────────────────────────────────────────────────────────
fig.text(0.5, 0.005,
         'Delta-delta interaction test (delta_CH vs delta_SN, Welch t-test, BH FDR p_adj < 0.05).  '
         'AW-M-3 excluded.  Chromatin: Keep + Review only.  Color = Log2 FC vs Naive.',
         ha='center', va='bottom',
         fontsize=7, fontfamily='Arial', color='#888888')

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = os.path.join(SCRIPT_DIR, 'Translocation_Summary_AllConditions.pdf')
out_png = os.path.join(SCRIPT_DIR, 'Translocation_Summary_AllConditions.png')

with PdfPages(out_pdf) as pdf:
    pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(out_png, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_png}')
