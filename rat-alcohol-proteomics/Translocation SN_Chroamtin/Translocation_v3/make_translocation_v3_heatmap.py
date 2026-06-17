"""
Translocation v3 — Heatmap matching existing panel style
=========================================================
Proteins : moved INTO Chromatin in AW (p_adj_BH < 0.05, from Translocation_v3_stats.xlsx)

Layout (matches make_heatmap_panel.py / make_heatmap_chrom_nuc_AW.py):
  Two stacked horizontal panels:
    Top    : Chromatin  — 3 thin rows (Intox, AW, PA)
    Bottom : Soluble Nuclear — 3 thin rows (Intox, AW, PA)
  Proteins as columns, sorted by Chromatin AW FC vs Naive (high → low)
  Color = Log2 Fold Change vs Naive
  Colormap: blue → white → red  (same as all other heatmaps)
  vmax = 2.0

Output
------
  Translocation_v3_heatmap.pdf
  Translocation_v3_heatmap.png
  Translocation_v3_heatmap_sourcedata.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STATS_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_v3_stats.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap.png')
OUT_CSV    = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_sourcedata.csv')

# ── Load significant proteins, sort by Chrom AW FC ───────────────────────────
sig = pd.read_excel(STATS_XLSX, sheet_name='Into_Chromatin')
sig = sig.sort_values('Chrom_FC_AW', ascending=False).reset_index(drop=True)
n_prot = len(sig)
print(f'Proteins moved into Chromatin (p_adj<0.05): {n_prot}')

if n_prot == 0:
    print('No significant proteins — no heatmap produced.')
    raise SystemExit

genes = sig['Gene'].tolist()

# FC matrices: shape (3, n_prot) — rows = Intox, AW, PA
CH_MAT = np.vstack([
    sig['Chrom_FC_Intox'].values,
    sig['Chrom_FC_AW'].values,
    sig['Chrom_FC_PA'].values,
])
SN_MAT = np.vstack([
    sig['SN_FC_Intox'].values,
    sig['SN_FC_AW'].values,
    sig['SN_FC_PA'].values,
])

COND_LABELS = ['Intox', 'AW', 'PA']

# ── Colormap (matches volcano plot / pub heatmap: blue → white → pink-red) ────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'],
    N=512
)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

COMPARTMENT_COLORS = {
    'Chromatin':       '#F7F3E3',
    'Soluble Nuclear': '#E3F2FA',
}

# ── Layout constants (same as make_heatmap_panel.py) ─────────────────────────
N_COMP   = 2
N_COND   = 3
ROW_H    = 0.55
BAND_GAP = 0.20
TITLE_H  = 0.40
TOP_PAD  = 0.55
BOT_PAD  = 1.10
LEFT_PAD = 2.20
RIGHT_PAD= 2.80

hm_w   = n_prot * 0.010          # inches — protein columns
FIG_W  = max(8.0, LEFT_PAD + hm_w + RIGHT_PAD)
BAND_H = N_COND * ROW_H
fig_h  = TOP_PAD + N_COMP * (TITLE_H + BAND_H) + (N_COMP - 1) * BAND_GAP + BOT_PAD

fig = plt.figure(figsize=(FIG_W, fig_h))
fig.patch.set_facecolor('white')

def ffy(y): return y / fig_h
def ffx(x): return x / FIG_W

panels = [
    ('Chromatin',       CH_MAT),
    ('Soluble Nuclear', SN_MAT),
]

im_ref = None

for idx, (label, mat) in enumerate(panels):
    band_top = fig_h - TOP_PAD - idx * (TITLE_H + BAND_H + BAND_GAP) - TITLE_H
    band_bot = band_top - BAND_H

    ax = fig.add_axes([ffx(LEFT_PAD), ffy(band_bot), ffx(hm_w), ffy(BAND_H)])
    ax.set_facecolor(COMPARTMENT_COLORS[label])

    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    if im_ref is None:
        im_ref = im

    # Row dividers
    for r in range(1, N_COND):
        ax.axhline(r - 0.5, color='white', linewidth=1.5, zorder=3)

    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Condition labels on right
    for r, cname in enumerate(COND_LABELS):
        y_frac = ffy(band_bot + BAND_H - (r + 0.5) * ROW_H)
        fig.text(ffx(LEFT_PAD + hm_w + 0.15), y_frac, cname,
                 ha='left', va='center', fontsize=10,
                 fontfamily='Arial', color='#111111')

    # Compartment title + n=
    fig.text(ffx(LEFT_PAD - 0.10), ffy(band_top + TITLE_H * 0.60),
             label, ha='left', va='center',
             fontsize=13, fontweight='bold', fontfamily='Arial', color='#111111')
    fig.text(ffx(LEFT_PAD - 0.10), ffy(band_top + TITLE_H * 0.12),
             f'n = {n_prot} proteins', ha='left', va='center',
             fontsize=8.5, fontfamily='Arial', color='#555555', style='italic')

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar_w = min(2.5, FIG_W * 0.25)
ax_cbar = fig.add_axes([ffx((FIG_W - cbar_w) / 2), ffy(0.28),
                         ffx(cbar_w), ffy(0.18)])
cb = fig.colorbar(im_ref, cax=ax_cbar, orientation='horizontal')
cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
cb.ax.tick_params(labelsize=8)
ax_cbar.text(0.5, -0.9, 'Log2 Fold Change vs Naive',
             transform=ax_cbar.transAxes,
             ha='center', va='top', fontsize=8.5,
             fontfamily='Arial', color='#444444')

fig.text(0.5, 0.01,
         'AW-M-3 excluded.  Sorted by Chromatin AW FC (high to low).  '
         'Interaction test: delta_CH vs delta_SN, Welch t-test BH FDR p_adj < 0.05',
         ha='center', va='bottom', fontsize=6.5,
         fontfamily='Arial', color='#666666')

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(0.5, 1 - (TOP_PAD * 0.15) / fig_h,
         f'Proteins moved into Chromatin during Acute Withdrawal  (n = {n_prot})',
         ha='center', va='top', fontsize=11, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')

# ── Save ──────────────────────────────────────────────────────────────────────
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')

# ── Source data ───────────────────────────────────────────────────────────────
src = pd.DataFrame({
    'Gene':            genes,
    'Chrom_FC_Intox':  CH_MAT[0],
    'Chrom_FC_AW':     CH_MAT[1],
    'Chrom_FC_PA':     CH_MAT[2],
    'SN_FC_Intox':     SN_MAT[0],
    'SN_FC_AW':        SN_MAT[1],
    'SN_FC_PA':        SN_MAT[2],
})
src.to_csv(OUT_CSV, index=False)
print(f'Saved: {OUT_CSV}')
