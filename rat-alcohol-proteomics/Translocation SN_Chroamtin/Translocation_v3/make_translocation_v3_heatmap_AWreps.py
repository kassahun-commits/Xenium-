"""
Translocation v3 — AW individual replicates heatmaps
=====================================================
Produces two PNG/PDF files:
  1. All 346 proteins (into chromatin, p_adj < 0.05)
  2. 60 true translocation proteins (CH up AND SN down, SN_FC_AW < -0.5)

Both show SN reps | Chromatin reps side by side, z-scored per protein.
Colormap: blue → white → pink-red (matches publication heatmaps).
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
FILE       = os.path.join(SCRIPT_DIR, '../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
STATS_XLSX = os.path.join(SCRIPT_DIR, 'Translocation_v3_stats.xlsx')

CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'],
    N=512
)
VMAX = 2.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

SN_COLS = ['Nuc_AW-F-1', 'Nuc_AW-F-2', 'Nuc_AW-M-1', 'Nuc_AW-M-2']
CH_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']

# ── Load raw data ──────────────────────────────────────────────────────────────
ch_raw = pd.read_excel(FILE, sheet_name='Chromatin')
ch_df  = ch_raw[ch_raw['Filter'].isin(['Keep', 'Review'])].copy()
ch_df['Gene symbol'] = ch_df['Gene symbol'].astype(str).str.strip()
ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')

sn_df = pd.read_excel(FILE, sheet_name='Soluble nuclear').copy()
sn_df['Gene symbol'] = sn_df['Gene symbol'].astype(str).str.strip()
sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')

# ── Load protein lists ─────────────────────────────────────────────────────────
all346 = pd.read_excel(STATS_XLSX, sheet_name='Into_Chromatin')
all346 = all346.sort_values('Interaction_score', ascending=False).reset_index(drop=True)

true60 = pd.read_excel(STATS_XLSX, sheet_name='G1_CH_up_SN_down_n60')
true60 = true60.sort_values('Interaction_score', ascending=False).reset_index(drop=True)

# ── Helper: build z-score matrix ──────────────────────────────────────────────
def build_matrix(genes):
    mat = np.full((len(genes), 8), np.nan)
    for i, gene in enumerate(genes):
        sn_vals = pd.to_numeric(sn_df.loc[gene, SN_COLS], errors='coerce').values if gene in sn_df.index else np.full(4, np.nan)
        ch_vals = pd.to_numeric(ch_df.loc[gene, CH_COLS], errors='coerce').values if gene in ch_df.index else np.full(4, np.nan)
        row = np.concatenate([sn_vals, ch_vals])
        mu  = np.nanmean(row)
        sd  = np.nanstd(row, ddof=1)
        mat[i] = (row - mu) / sd if sd > 0 else row - mu
    return mat

# ── Helper: draw and save heatmap ─────────────────────────────────────────────
def make_heatmap(genes, mat, title, subtitle, out_pdf, out_png):
    n_prot   = len(genes)
    ROW_H_IN = 0.018
    fig_h    = max(6.0, n_prot * ROW_H_IN + 2.5)
    FIG_W    = max(5.0, n_prot * 0.06 + 1.5)

    fig, ax = plt.subplots(figsize=(FIG_W, fig_h))
    fig.patch.set_facecolor('white')

    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax.axvline(3.5, color='white', linewidth=2.5, zorder=3)

    col_labels = ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2',
                  'AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2']
    ax.set_xticks(range(8))
    ax.set_xticklabels(col_labels, fontsize=7.5, fontfamily='Arial',
                       rotation=45, ha='left')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='x', which='both', length=0, pad=2)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.annotate('Soluble Nuclear', xy=(1.5/8, 1.0), xycoords='axes fraction',
                xytext=(1.5/8, 1.18), textcoords='axes fraction',
                ha='center', va='bottom', fontsize=9, fontweight='bold',
                fontfamily='Arial', annotation_clip=False)
    ax.annotate('Chromatin', xy=(5.5/8, 1.0), xycoords='axes fraction',
                xytext=(5.5/8, 1.18), textcoords='axes fraction',
                ha='center', va='bottom', fontsize=9, fontweight='bold',
                fontfamily='Arial', annotation_clip=False)

    cbar = fig.colorbar(im, ax=ax, orientation='horizontal',
                        pad=0.12, fraction=0.03, shrink=0.6)
    cbar.set_ticks([-2, -1, 0, 1, 2])
    cbar.ax.tick_params(labelsize=7.5)
    cbar.set_label('Z-score (LFQ intensity)', fontsize=8.5, fontfamily='Arial')

    fig.text(0.5, 0.005,
             'AW-M-3 excluded.  Sorted by interaction score (high to low).',
             ha='center', va='bottom', fontsize=6.5,
             fontfamily='Arial', color='#666666')

    ax.set_title(f'{title}  (n = {n_prot} proteins)\n{subtitle}',
                 fontsize=9.5, fontweight='bold', fontfamily='Arial',
                 color='#1A3A5C', pad=95)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_png}')

# ── Heatmap 1: all 346 ────────────────────────────────────────────────────────
genes_346 = all346['Gene'].tolist()
mat_346   = build_matrix(genes_346)
make_heatmap(
    genes_346, mat_346,
    title    = 'Into Chromatin during AW — all proteins',
    subtitle = 'CH increase > SN increase vs naive, BH FDR p_adj < 0.05',
    out_pdf  = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_AWreps_346.pdf'),
    out_png  = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_AWreps_346.png'),
)

# ── Heatmap 2: 60 true translocation ─────────────────────────────────────────
genes_60 = true60['Gene'].tolist()
mat_60   = build_matrix(genes_60)
make_heatmap(
    genes_60, mat_60,
    title    = 'True translocation into Chromatin during AW',
    subtitle = 'CH up AND SN down vs naive, BH FDR p_adj < 0.05',
    out_pdf  = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_AWreps_60.pdf'),
    out_png  = os.path.join(SCRIPT_DIR, 'Translocation_v3_heatmap_AWreps_60.png'),
)
