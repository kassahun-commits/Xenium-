"""
Chrom/SN Balance Ratio Heatmap — all 4 conditions
===================================================
Uses the pre-computed Chrom_vs_SN_Balance.xlsx.

For each protein, shows the actual Chrom/SN ratio for each condition:
  Naive  = RATIO-Naive
  Intox  = RATIO-Naive + FC_Intox
  AW     = RATIO-Naive + FC_AW
  PA     = RATIO-Naive + FC_PA

Significance already computed in the file (Sig column = Up/Down/NS,
corrected p-value > 3.3 threshold).

Shows proteins significant in ≥ 1 condition (Intox / AW / PA vs Naive).
Proteins sorted by AW ratio (high → low).

Same VP colormap + style as Translocation_v3.

Outputs (saved in same folder as this script):
  Balance_Ratio_sig_proteins.xlsx
  Balance_Ratio_heatmap.pdf / .png
"""

import os
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_XLSX    = os.path.join(SCRIPT_DIR, 'Chrom_vs_SN_Balance.xlsx')
OUT_XLSX   = os.path.join(SCRIPT_DIR, 'Balance_Ratio_sig_proteins.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'Balance_Ratio_heatmap.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Balance_Ratio_heatmap.png')

CORR_THRESH = 3.3

# ── Load ───────────────────────────────────────────────────────────────────────
print('Loading data...')
df = pd.read_excel(IN_XLSX, sheet_name='Chrom_vs_SN_Balance')
df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

# ── Compute actual ratio per condition ─────────────────────────────────────────
df['Ratio_Naive'] = pd.to_numeric(df['RATIO-Naive'], errors='coerce')
df['Ratio_Intox'] = df['Ratio_Naive'] + pd.to_numeric(df['FC_Intox'], errors='coerce')
df['Ratio_AW']    = df['Ratio_Naive'] + pd.to_numeric(df['FC_AW'],    errors='coerce')
df['Ratio_PA']    = df['Ratio_Naive'] + pd.to_numeric(df['FC_PA'],    errors='coerce')

# ── Filter: significant in ≥ 1 condition ──────────────────────────────────────
sig_mask = (df['Sig_Intox'] != 'NS') | (df['Sig_AW'] != 'NS') | (df['Sig_PA'] != 'NS')
sig      = df[sig_mask].copy()
print(f'Proteins sig in ≥1 condition: {len(sig)}')
print(f'  Intox Up/Down: {(sig["Sig_Intox"]=="Up").sum()}/{(sig["Sig_Intox"]=="Down").sum()}')
print(f'  AW    Up/Down: {(sig["Sig_AW"]=="Up").sum()}/{(sig["Sig_AW"]=="Down").sum()}')
print(f'  PA    Up/Down: {(sig["Sig_PA"]=="Up").sum()}/{(sig["Sig_PA"]=="Down").sum()}')

# Sort by AW ratio (high → low) for display
sig = sig.sort_values('Ratio_AW', ascending=False, na_position='last').reset_index(drop=True)

# ── Save Excel ─────────────────────────────────────────────────────────────────
out_cols = ['Accession', 'Gene symbol', 'Description', 'Filter',
            'Ratio_Naive', 'Ratio_Intox', 'Ratio_AW', 'Ratio_PA',
            'FC_Intox', 'Corrected_Intox', 'Sig_Intox',
            'FC_AW',    'Corrected_AW',    'Sig_AW',
            'FC_PA',    'Corrected_PA',    'Sig_PA']
sig[out_cols].to_excel(OUT_XLSX, index=False)
print(f'Saved: {OUT_XLSX}')

# ── Build heatmap matrix ───────────────────────────────────────────────────────
CONDS       = ['Naive', 'Intox', 'AW', 'PA']
COND_LABELS = ['Naïve', 'Intoxication', 'Acute\nWithdrawal', 'Protracted\nAbstinence']
RATIO_COLS  = ['Ratio_Naive', 'Ratio_Intox', 'Ratio_AW', 'Ratio_PA']
SIG_COLS    = ['Sig_Naive',   'Sig_Intox',   'Sig_AW',   'Sig_PA']

# Sig_Naive not really used (it's the baseline) — leave blank
sig['Sig_Naive'] = 'NS'

n_prots = len(sig)
n_conds = 4

mat = np.full((n_prots, n_conds), np.nan)
for ci, col in enumerate(RATIO_COLS):
    mat[:, ci] = pd.to_numeric(sig[col], errors='coerce').values

# Significance annotation strip (per condition, skip Naive)
# 1 = Up (red), 0 = NS (gray), -1 = Down (blue)
sig_mat = np.zeros((n_prots, n_conds), dtype=float)
for ci, col in enumerate(SIG_COLS):
    vals = sig[col].fillna('NS')
    sig_mat[:, ci] = np.where(vals == 'Up', 1, np.where(vals == 'Down', -1, 0))

# ── Figure sizing (same logic as v3) ──────────────────────────────────────────
INCHES_PER_PROT = 0.006
ROW_H           = max(0.006, min(0.025, INCHES_PER_PROT))
ANNOT_ROW_H     = 0.10
TOP_MARGIN      = 1.4
BOT_MARGIN      = 0.7
LEFT_MARGIN     = 0.8
RIGHT_MARGIN    = 2.0
COL_W           = 0.55    # inches per condition column

HEAT_H   = n_prots * ROW_H
HEAT_W   = n_conds * COL_W
FIG_H    = TOP_MARGIN + ANNOT_ROW_H * 3 + HEAT_H + BOT_MARGIN
FIG_W    = LEFT_MARGIN + HEAT_W + RIGHT_MARGIN

print(f'Figure size: {FIG_W:.1f} × {FIG_H:.1f} inches  ({n_prots} proteins)')

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

# Axes positions in figure fraction
def frac(inches, total): return inches / total

ax_l = frac(LEFT_MARGIN, FIG_W)
ax_r = frac(LEFT_MARGIN + HEAT_W, FIG_W)
ax_b = frac(BOT_MARGIN, FIG_H)
ax_t = frac(BOT_MARGIN + HEAT_H, FIG_H)

# ── Annotation strips (3 rows: Intox, AW, PA) above heatmap ───────────────────
ANNOT_CMAP = ListedColormap(['#2B7FD4', '#E8E8E8', '#E8305A'])
ANNOT_NORM = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ANNOT_CMAP.N)

annot_axes = []
for ai in range(3):   # Intox=0, AW=1, PA=2
    a_b = frac(BOT_MARGIN + HEAT_H + ai * ANNOT_ROW_H, FIG_H)
    a_t = frac(BOT_MARGIN + HEAT_H + (ai + 1) * ANNOT_ROW_H, FIG_H)
    aax = fig.add_axes([ax_l, a_b, ax_r - ax_l, a_t - a_b])
    aax.imshow(sig_mat[:, ai + 1].reshape(1, -1),  # skip Naive col
               aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM,
               interpolation='nearest')
    aax.set_xticks([])
    aax.set_yticks([])
    for sp in aax.spines.values(): sp.set_visible(False)
    annot_axes.append(aax)

# Label the annotation strip rows
strip_labels = ['Intox', 'AW', 'PA']
for ai, (aax, lbl) in enumerate(zip(annot_axes, strip_labels)):
    aax.text(-0.015, 0.5, lbl,
             transform=aax.transAxes, ha='right', va='center',
             fontsize=7, fontfamily='Arial', color='#444444')

# ── Main heatmap ───────────────────────────────────────────────────────────────
ax = fig.add_axes([ax_l, ax_b, ax_r - ax_l, ax_t - ax_b])

VP_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'VP', ['#1A5FA0', '#FFFFFF', '#C01E42'])

# Set colormap limits based on data
abs_max = np.nanpercentile(np.abs(mat), 98)
vmax    = min(abs_max, 8)
norm    = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

im = ax.imshow(mat, aspect='auto', cmap=VP_CMAP, norm=norm,
               interpolation='nearest')

# Vertical lines between conditions
for xi in range(1, n_conds):
    ax.axvline(xi - 0.5, color='#AAAAAA', lw=0.5, zorder=2)

ax.set_xticks(range(n_conds))
ax.set_xticklabels(COND_LABELS, fontsize=10, fontfamily='Arial')
ax.set_yticks([])
ax.tick_params(axis='x', length=0, pad=4)
for sp in ax.spines.values(): sp.set_visible(False)

# n label on left
ax.text(-0.015, 0.5, f'n = {n_prots}',
        transform=ax.transAxes, ha='right', va='center',
        fontsize=9, fontfamily='Arial', color='#444444',
        rotation=90)

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar_l = ax_r + frac(0.12, FIG_W)
cbar_b = ax_b + (ax_t - ax_b) * 0.25
cbar_w = frac(0.12, FIG_W)
cbar_h = (ax_t - ax_b) * 0.50

cax = fig.add_axes([cbar_l, cbar_b, cbar_w, cbar_h])
cb  = fig.colorbar(im, cax=cax, extend='both')
cb.set_label('Chrom/SN Ratio\n(Log2)', fontsize=9, fontfamily='Arial', labelpad=6)
cb.ax.tick_params(labelsize=8)
ticks = [-vmax, -vmax/2, 0, vmax/2, vmax]
cb.set_ticks([t for t in ticks if abs(t) <= vmax])

# ── Annotation legend ──────────────────────────────────────────────────────────
leg_x = cbar_l
leg_y = ax_b + (ax_t - ax_b) * 0.10
handles = [
    mpatches.Patch(color='#E8305A', label='Chrom enriched (Up)'),
    mpatches.Patch(color='#2B7FD4', label='SN enriched (Down)'),
    mpatches.Patch(color='#E8E8E8', label='NS'),
]
fig.legend(handles=handles, loc='lower right',
           bbox_to_anchor=(0.99, 0.02),
           fontsize=8, framealpha=0.9,
           edgecolor='#CCCCCC',
           title='Significance strip\n(vs Naïve)', title_fontsize=8)

# ── Title ─────────────────────────────────────────────────────────────────────
top_y = frac(BOT_MARGIN + HEAT_H + 3 * ANNOT_ROW_H + 0.15, FIG_H)
fig.text(0.5, top_y + frac(0.4, FIG_H),
         f'Chrom/SN Balance Ratio — Proteins Significant in ≥1 Condition  (n = {n_prots})',
         ha='center', va='bottom', fontsize=13, fontweight='bold', fontfamily='Arial')

fig.text(0.5, frac(0.1, FIG_H),
         f'Sorted by AW ratio (high→low)  |  Corrected p > {CORR_THRESH}  |  '
         'Positive = Chrom-enriched; Negative = SN-enriched  |  Keep+Review filter applied',
         ha='center', va='bottom', fontsize=8, color='#666666',
         style='italic', fontfamily='Arial')

# ── Save ──────────────────────────────────────────────────────────────────────
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
