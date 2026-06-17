"""
Translocation v4 — Chromatin / Soluble Nuclear log2 ratio heatmap
=================================================================
For every protein present in both fractions (Chromatin Keep+Review, SN all)
compute the log2(Chrom/SN) ratio per condition:

    Ratio = mean(log2 Chrom reps) − mean(log2 SN reps)

Four conditions: Naive, Intoxication, Acute Withdrawal (AW-M-3 excluded), PA

Heatmap:  proteins as columns, 4 condition rows, single band
          VP diverging colormap (same as Translocation v3)
          Sorted by Naive ratio high → low

Outputs:
  Translocation_v4_stats.xlsx
  Translocation_v4_heatmap.pdf/.png
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
FILE = os.path.join(SCRIPT_DIR,
       '../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_v4_stats.xlsx')

# ── Column groups (AW-M-3 excluded from AW) ──────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',   'Nuc_N-F-2',   'Nuc_N-F-3',   'Nuc_N-M-1',   'Nuc_N-M-2']

CH_I_COLS  = ['Chrom_I-F-1', 'Chrom_I-F-2', 'Chrom_I-F-3', 'Chrom_I-M-1', 'Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',   'Nuc_I-F-2',   'Nuc_I-F-3',   'Nuc_I-M-1',   'Nuc_I-M-2']

CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']   # M-3 excluded
SN_AW_COLS = ['Nuc_AW-F-1',   'Nuc_AW-F-2',   'Nuc_AW-M-1',   'Nuc_AW-M-2']

CH_PA_COLS = ['Chrom_PA-F-1', 'Chrom_PA-F-2', 'Chrom_PA-F-3', 'Chrom_PA-M-1', 'Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1',   'Nuc_PA-F-2',   'Nuc_PA-F-3',   'Nuc_PA-M-1',   'Nuc_PA-M-2']

CONDITIONS = [
    ('Naive',                CH_N_COLS,  SN_N_COLS),
    ('Intoxication',         CH_I_COLS,  SN_I_COLS),
    ('Acute Withdrawal',     CH_AW_COLS, SN_AW_COLS),
    ('Protracted Abstinence',CH_PA_COLS, SN_PA_COLS),
]

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

# ── Heatmap layout ────────────────────────────────────────────────────────────
ROW_H           = 0.35      # inches per condition row
N_COND          = 4
BAND_H          = N_COND * ROW_H
LEFT_IN         = 2.20
RIGHT_IN        = 2.80
TOP_PAD         = 0.55
BOT_PAD         = 1.20
CBAR_W_IN       = 2.8
CBAR_H_IN       = 0.14
CBAR_BOT        = 0.48
BAND_BG         = '#FBF9F0'   # same warm tint as v3 Chromatin band

COND_LABELS = ['Naive', 'Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

# ── Helper ────────────────────────────────────────────────────────────────────
def col_mean(df, gene, cols):
    vals = pd.to_numeric(df.loc[gene, cols], errors='coerce').values \
           if gene in df.index else np.full(len(cols), np.nan)
    return float(np.nanmean(vals))

# ── Load data ─────────────────────────────────────────────────────────────────
print('Loading Chromatin sheet (Keep + Review)...')
ch_raw = pd.read_excel(FILE, sheet_name='Chromatin')
ch_df  = ch_raw[ch_raw['Filter'].isin(['Keep', 'Review'])].copy().reset_index(drop=True)

print('Loading Soluble Nuclear sheet...')
sn_df  = pd.read_excel(FILE, sheet_name='Soluble nuclear').copy().reset_index(drop=True)

for df in [ch_df, sn_df]:
    df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')

common_genes = sorted(set(ch_df.index) & set(sn_df.index))
print(f'Proteins in both fractions: {len(common_genes)}')

# ── Compute log2(Chrom/SN) ratios ─────────────────────────────────────────────
rows = []
for gene in common_genes:
    row = {'Gene': gene}
    for cond_name, ch_cols, sn_cols in CONDITIONS:
        ch_mean = col_mean(ch_df, gene, ch_cols)
        sn_mean = col_mean(sn_df, gene, sn_cols)
        # log2(Chrom/SN) = log2_Chrom - log2_SN (values already in log2 space)
        ratio = (ch_mean - sn_mean) if (np.isfinite(ch_mean) and np.isfinite(sn_mean)) else np.nan
        col_label = f'Ratio_{cond_name.replace(" ", "_")}'
        row[col_label] = round(ratio, 4) if np.isfinite(ratio) else np.nan
    rows.append(row)

df_out = pd.DataFrame(rows)

# Sort by Naive ratio high → low (proteins most enriched in Chromatin first)
df_out = df_out.sort_values('Ratio_Naive', ascending=False).reset_index(drop=True)

print(f'Proteins with valid Naive ratio: {df_out["Ratio_Naive"].notna().sum()}')
print(f'Ratio_Naive  range: {df_out["Ratio_Naive"].min():.2f} to {df_out["Ratio_Naive"].max():.2f}')
print(f'Ratio_Intoxication range: {df_out["Ratio_Intoxication"].min():.2f} to {df_out["Ratio_Intoxication"].max():.2f}')

# ── Save Excel ─────────────────────────────────────────────────────────────────
with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    df_out.to_excel(writer, sheet_name='Chrom_SN_Ratios', index=False)
print(f'Saved: {OUT_XLSX}')

# ── Heatmap ───────────────────────────────────────────────────────────────────
ratio_cols = [f'Ratio_{c.replace(" ", "_")}' for _, c, _ in
              [(n, n, None) for n, _, _ in CONDITIONS]]

# Build matrix: shape (4, n_proteins), drop rows where all NaN
df_hm = df_out.dropna(subset=ratio_cols, how='all').reset_index(drop=True)
n_prot = len(df_hm)
print(f'Proteins in heatmap: {n_prot}')

mat = np.vstack([df_hm[col].fillna(0).values for col in ratio_cols])

inches_per_prot = 0.004
hm_w   = max(n_prot * inches_per_prot, 3.0)
FIG_W  = LEFT_IN + hm_w + RIGHT_IN
FIG_H  = TOP_PAD + BAND_H + BOT_PAD

fig      = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
y_bot    = BOT_PAD
y_top    = BOT_PAD + BAND_H

ax = fig.add_axes([LEFT_IN / FIG_W, y_bot / FIG_H,
                   hm_w / FIG_W,    BAND_H / FIG_H],
                  facecolor=BAND_BG)

im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')

# Row dividers between conditions
for r in range(1, N_COND):
    ax.axhline(r - 0.5, color='white', linewidth=1.2, zorder=3)

ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)

# Condition labels on the right
for r, cname in enumerate(COND_LABELS):
    y_frac = (y_bot + BAND_H - (r + 0.5) * ROW_H) / FIG_H
    fig.text((LEFT_IN + hm_w + 0.12) / FIG_W, y_frac,
             cname, ha='left', va='center',
             fontsize=8, fontfamily='Arial', color='#222222')

# Band label on the left
mid_y = (y_bot + BAND_H / 2) / FIG_H
fig.text((LEFT_IN - 0.12) / FIG_W, mid_y,
         'Chromatin / SN',
         ha='right', va='center',
         fontsize=10, fontweight='bold',
         fontfamily='Arial', color='#111111')
fig.text((LEFT_IN - 0.12) / FIG_W, mid_y - 0.06,
         f'log2 ratio',
         ha='right', va='top',
         fontsize=7.5, fontfamily='Arial',
         color='#555555', style='italic')
fig.text((LEFT_IN - 0.12) / FIG_W, mid_y - 0.10,
         f'n = {n_prot}',
         ha='right', va='top',
         fontsize=7.5, fontfamily='Arial',
         color='#555555', style='italic')

# Colorbar
cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H,
                       CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
cb.ax.tick_params(labelsize=7)
ax_cb.set_xlabel('log2(Chromatin / Soluble Nuclear)',
                 fontsize=8, fontfamily='Arial', color='#444444', labelpad=4)

# Title
fig.text(0.5, (FIG_H - TOP_PAD * 0.35) / FIG_H,
         f'Chromatin / Soluble Nuclear Ratio — All Conditions  (n = {n_prot})',
         ha='center', va='top',
         fontsize=9, fontweight='bold',
         fontfamily='Arial', color='#1A3A5C')

# Footnote
fig.text(0.5, 0.01,
         'log2(Chrom/SN) = mean(log2 Chrom reps) − mean(log2 SN reps).  '
         'Chromatin: Keep + Review only.  AW-M-3 excluded.  '
         'Sorted by Naive ratio high → low.',
         ha='center', va='bottom', fontsize=5.5,
         fontfamily='Arial', color='#888888')

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = os.path.join(SCRIPT_DIR, 'Translocation_v4_heatmap.pdf')
out_png = os.path.join(SCRIPT_DIR, 'Translocation_v4_heatmap.png')

with PdfPages(out_pdf) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out_png}')
