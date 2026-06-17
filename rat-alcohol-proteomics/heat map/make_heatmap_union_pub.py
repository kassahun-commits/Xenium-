"""
Publication-quality union heatmap panel — 4 compartments.
- Uses EDIT Excluding AWM3 file (AW-M-3 already excluded in data)
- Chromatin: Keep + Review filter only; all others: unfiltered
- Union proteins: corrected p-value >= 3.3 in ANY of the 3 conditions
- 3 condition rows per compartment: Intoxication, Acute Withdrawal, Protracted Abstinence
- Fold-change columns used directly from the precomputed spreadsheet
- Width of each compartment band proportional to protein count
- Colormap matches volcano plot colors: blue (#2B7FD4) → white → pink-red (#E8305A)
- Proteins within each band sorted by Acute Withdrawal fold change
"""

import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patheffects as pe

# ── File & thresholds ──────────────────────────────────────────────────────────
FILE       = '../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
CORR_THRESH = 3.3
VMAX        = 3.0        # symmetric color-scale clip (Log2 FC)

# ── Colormap: VP blue → white → VP pink-red ───────────────────────────────────
# Matches volcano plot colors: C_DOWN=#2B7FD4, C_UP=#E8305A
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    [
        '#1A5FA0',   # deep blue
        '#2B7FD4',   # VP blue (C_DOWN)
        '#6AAEE0',   # mid blue
        '#AECFE8',   # light blue
        '#DDEEF8',   # very light blue
        '#FFFFFF',   # white (centre)
        '#FDD8E7',   # very light pink
        '#F5A0BC',   # light pink
        '#EE5F8B',   # mid pink
        '#E8305A',   # VP pink-red (C_UP)
        '#C01E42',   # deep red
    ],
    N=512,
)

# ── Compartments: (display label, sheet name, filter mode) ────────────────────
COMPARTMENTS = [
    ('Membrane',        'Membrane',        'all'),
    ('Cytosol',         'Cytosol',         'all'),
    ('Chromatin',       'Chromatin',       'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'all'),
]

# Column mapping: (FC col, Corrected col, condition label)
CONDITIONS = [
    ('Fold change',   'Corrected',   'Intoxication'),
    ('Fold change.1', 'Corrected.1', 'Acute Withdrawal'),
    ('Fold change.2', 'Corrected.2', 'Protracted Abstinence'),
]

# ── Load & filter data ─────────────────────────────────────────────────────────
panel_data = []
for disp, sheet, fmode in COMPARTMENTS:
    df = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in df.columns:
        df = df[df['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    # Union: significant in at least one condition
    union_mask = pd.Series(False, index=df.index)
    for fc_col, corr_col, _ in CONDITIONS:
        sig = df[corr_col].notna() & (df[corr_col] >= CORR_THRESH)
        union_mask = union_mask | sig

    df_u = df[union_mask].reset_index(drop=True)

    # Build (n_conditions × n_proteins) fold-change matrix
    rows = []
    for fc_col, _, _ in CONDITIONS:
        rows.append(pd.to_numeric(df_u[fc_col], errors='coerce').fillna(0).values)
    mat = np.vstack(rows)   # shape (3, n_proteins)

    # Sort by AW fold change (row index 1)
    order = np.argsort(mat[1])
    mat   = mat[:, order]

    n = len(df_u)
    panel_data.append(dict(label=disp, n=n, mat=mat))
    print(f'{disp}: {n} proteins in union')

# ── Layout constants ───────────────────────────────────────────────────────────
INCHES_PER_PROT  = 0.0072   # width per protein column
ROW_H            = 0.50     # height per condition row (inches)
N_COND           = len(CONDITIONS)
BAND_H           = N_COND * ROW_H   # total band height per compartment

LEFT_IN          = 2.40     # space left of heatmap (for compartment label)
RIGHT_IN         = 2.70     # space right of heatmap (for condition labels)
TOP_PAD          = 0.55     # above first band
BOT_PAD          = 0.90     # below last band (colorbar area)
SECTION_GAP      = 0.32     # vertical gap between bands
LABEL_OFFSET     = 0.14     # spacing between label and heatmap edge

CBAR_W_IN        = 2.8
CBAR_H_IN        = 0.14
CBAR_BOT         = 0.22

# ── Figure geometry ────────────────────────────────────────────────────────────
max_band_w = max(d['n'] * INCHES_PER_PROT for d in panel_data)
FIG_W = LEFT_IN + max_band_w + RIGHT_IN
FIG_H = (TOP_PAD
         + len(panel_data) * BAND_H
         + (len(panel_data) - 1) * SECTION_GAP
         + BOT_PAD)

# ── Condition label short forms ────────────────────────────────────────────────
COND_SHORT = {
    'Intoxication':          'Intoxication',
    'Acute Withdrawal':      'Acute Withdrawal',
    'Protracted Abstinence': 'Protracted Abstinence',
}

# Subtle alternating background tones per compartment
COMP_BG = {
    'Membrane':        '#F7F3FA',
    'Cytosol':         '#F2F9F4',
    'Chromatin':       '#FBF9F0',
    'Soluble Nuclear': '#F0F6FB',
}

norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

im_ref = None
y_cursor = FIG_H - TOP_PAD   # starts from top of figure, moves down

for idx, d in enumerate(panel_data):
    label  = d['label']
    n_prot = d['n']
    mat    = d['mat']
    hm_w   = n_prot * INCHES_PER_PROT

    y_band_top = y_cursor
    y_band_bot = y_band_top - BAND_H

    # Figure-fraction coordinates for the heatmap axes
    ax_left   = LEFT_IN / FIG_W
    ax_bottom = y_band_bot / FIG_H
    ax_w      = hm_w / FIG_W
    ax_h      = BAND_H / FIG_H

    # ── Subtle section background (extends slightly beyond heatmap) ───────────
    bg_pad  = 0.06    # inches padding around the heatmap
    bg_rect = plt.matplotlib.patches.FancyBboxPatch(
        ((LEFT_IN - bg_pad) / FIG_W, (y_band_bot - bg_pad) / FIG_H),
        (hm_w + 2 * bg_pad) / FIG_W,
        (BAND_H + 2 * bg_pad) / FIG_H,
        boxstyle='round,pad=0.002',
        facecolor=COMP_BG.get(label, '#F8F8F8'),
        edgecolor='none',
        transform=fig.transFigure,
        zorder=0,
        clip_on=False,
    )
    fig.add_artist(bg_rect)

    # ── Heatmap axes ──────────────────────────────────────────────────────────
    ax = fig.add_axes([ax_left, ax_bottom, ax_w, ax_h])
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                   interpolation='nearest', rasterized=True)
    im_ref = im

    # Thin white horizontal lines between condition rows
    for r in range(1, N_COND):
        ax.axhline(r - 0.5, color='white', linewidth=1.6, zorder=3)

    # Minimal border
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.6)
        sp.set_color('#BBBBBB')
    ax.set_xticks([])
    ax.set_yticks([])

    # ── Compartment label — left of band, vertically centred ─────────────────
    cx_label = (LEFT_IN - LABEL_OFFSET) / FIG_W
    cy_label = (y_band_bot + BAND_H / 2) / FIG_H

    # Bold compartment name + protein count on next line
    fig.text(
        cx_label, cy_label,
        label,
        ha='right', va='bottom',
        fontsize=11.5, fontweight='bold',
        fontfamily='Arial', color='#1A1A1A',
    )
    fig.text(
        cx_label, cy_label,
        f'n = {n_prot}',
        ha='right', va='top',
        fontsize=8.5,
        fontfamily='Arial', color='#666666',
    )

    # ── Condition labels — right of band, one per row ─────────────────────────
    for r, (_, _, cname) in enumerate(CONDITIONS):
        y_row_center = (y_band_top - (r + 0.5) * ROW_H) / FIG_H
        fig.text(
            (LEFT_IN + hm_w + LABEL_OFFSET) / FIG_W,
            y_row_center,
            COND_SHORT[cname],
            ha='left', va='center',
            fontsize=9, fontfamily='Arial', color='#333333',
        )

    y_cursor = y_band_bot - SECTION_GAP

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar_l = (FIG_W - CBAR_W_IN) / 2
ax_cbar = fig.add_axes([
    cbar_l / FIG_W,
    CBAR_BOT / FIG_H,
    CBAR_W_IN / FIG_W,
    CBAR_H_IN / FIG_H,
])
cb = fig.colorbar(im_ref, cax=ax_cbar, orientation='horizontal', extend='both')
cb.ax.tick_params(labelsize=8, length=3, width=0.6, color='#666666', labelcolor='#333333')
cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
cb.outline.set_linewidth(0.6)
cb.outline.set_edgecolor('#BBBBBB')

# Colorbar title above bar
ax_cbar.set_title(
    'Log₂ Fold Change vs Naive',
    fontsize=8.5, fontfamily='DejaVu Sans', color='#444444',
    pad=5,
)

# ── Footer note ───────────────────────────────────────────────────────────────
fig.text(
    0.5, 0.005,
    'Union: corrected p (-log2) > 3.3 in at least one condition  |  '
    'Chromatin: Keep & Review only  |  AW-M-3 excluded',
    ha='center', va='bottom',
    fontsize=7.5, fontfamily='Arial', color='#888888', style='italic',
)

# ── Save ──────────────────────────────────────────────────────────────────────
OUT = 'Heatmap_Union_PubQuality.pdf'
pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
pdf.close()
plt.close(fig)
print(f'\nSaved: {OUT}')
print(f'Figure size: {FIG_W:.1f} × {FIG_H:.1f} inches')
