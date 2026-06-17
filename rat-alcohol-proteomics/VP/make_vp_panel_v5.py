"""
VP panel figure — style matches individual VPs exactly.
Uses the same dot sizes, fonts, and axes as make_vp_all12_final.py
but shows BOTH increased (red) AND decreased (blue) in every panel.
No jitter. AW-M-3 excluded. Chromatin: Keep+Review filter.

Output: VP_Panel_Figure_v5_combined.pdf
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'VP_Panel_Figure_v5_combined.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
XLIM = (-9.00, 9.74)
YLIM = (-0.5,  27.67)

C_UP   = '#E8305A'
C_DOWN = '#2B7FD4'
C_NS   = '#EBEBEB'

ROW_COLORS = {
    'Membrane':        '#E8E4F0',
    'Cytosol':         '#D6EDDC',
    'Chromatin':       '#F5F0D4',
    'Soluble Nuclear': '#D4E8F5',
}

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'Memb_N-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'Cyto_N-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'Chrom_N-', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'Nuc_N-',   'all'),
]

COND_SUFFIXES = {
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2'],
    'Protracted Abstinence': ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
CONDITIONS    = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
COL_TITLES    = CONDITIONS

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pv    = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pv.append(stats.ttest_ind(n, c, equal_var=False)[1]
                  if len(n) >= 2 and len(c) >= 2 else np.nan)
    pv    = pd.Series(pv, index=df.index)
    valid = pv.notna()
    ranks = pv[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pv[valid] * valid.sum() / ranks))
    return fc, corr

def draw_panel(ax, fc, corrected):
    """Draw one VP panel — identical style to individual VPs, both colors shown."""
    mask = fc.notna() & corrected.notna()
    f    = fc[mask].values
    cr   = corrected[mask].values

    up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
    down = (cr > CORR_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down
    n_up   = int(up.sum())
    n_down = int(down.sum())

    # Dot sizes + alphas match the individual VP exactly
    ax.scatter(f[ns],   cr[ns],   s=10, color=C_NS,   alpha=0.40, linewidths=0, rasterized=True)
    ax.scatter(f[down], cr[down], s=35, color=C_DOWN, alpha=0.90, linewidths=0, rasterized=True)
    ax.scatter(f[up],   cr[up],   s=35, color=C_UP,   alpha=0.90, linewidths=0, rasterized=True)

    # Threshold lines — same style as individual VP
    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlim(XLIM); ax.set_ylim(YLIM)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    # Font sizes match individual VP
    ax.set_xlabel('Log$_2$ Fold Change',          fontsize=16, labelpad=6)
    ax.set_ylabel('Corrected p-value\n(−log$_2$)', fontsize=16, labelpad=6)
    ax.tick_params(labelsize=12)

    ax.legend(handles=[
        mpatches.Patch(color=C_UP,   label=f'Increased  (n={n_up})'),
        mpatches.Patch(color=C_DOWN, label=f'Decreased  (n={n_down})'),
        mpatches.Patch(color=C_NS,   label='NS'),
    ], fontsize=11, loc='upper left', framealpha=0.88,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.7,
       handlelength=1.0, labelspacing=0.4)

# ── Pre-load data ──────────────────────────────────────────────────────────────
all_data = {}
for ri, (disp_label, sheet, prefix, n_pat, fmode) in enumerate(SHEETS):
    raw = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)
    naive_cols = [c for c in df.columns if n_pat in str(c)]
    for ci, cond_name in enumerate(CONDITIONS):
        suffixes  = COND_SUFFIXES[cond_name]
        cond_cols = [c for c in df.columns if any(prefix+s in str(c) for s in suffixes)]
        fc, corrected = calc_stats(df, naive_cols, cond_cols)
        all_data[(ri, ci)] = (fc, corrected, disp_label)
    print(f'Loaded {disp_label}')

# ── Figure — large enough that each panel matches individual VP quality ────────
N_ROWS, N_COLS = 4, 3

# Each panel should be ~8×6 inches like the individual VPs
# 3 cols × 8 = 24 + margins → 32 wide
# 4 rows × 6 = 24 + margins → 30 tall
fig = plt.figure(figsize=(32, 30))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(
    N_ROWS, N_COLS,
    figure  = fig,
    hspace  = 0.40,
    wspace  = 0.30,
    top     = 0.93,
    bottom  = 0.05,
    left    = 0.22,
    right   = 0.97,
)

axes_grid = {}
for ri in range(N_ROWS):
    for ci in range(N_COLS):
        fc, corrected, disp_label = all_data[(ri, ci)]
        ax = fig.add_subplot(gs[ri, ci])
        axes_grid[(ri, ci)] = ax
        draw_panel(ax, fc, corrected)

# ── Column headers ─────────────────────────────────────────────────────────────
for ci, title in enumerate(COL_TITLES):
    pos = axes_grid[(0, ci)].get_position()
    fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.015, title,
             ha='center', va='bottom',
             fontsize=22, fontweight='bold', fontfamily='Arial', color='#111111')

# ── Row labels — colored strip ─────────────────────────────────────────────────
for ri, (disp_label, *_) in enumerate(SHEETS):
    pos0 = axes_grid[(ri, 0)].get_position()

    box_right = pos0.x0 - 0.075
    box_left  = box_right - 0.040
    box_bot   = pos0.y0
    box_top   = pos0.y1

    box = FancyBboxPatch(
        (box_left, box_bot), box_right - box_left, box_top - box_bot,
        boxstyle='round,pad=0.002', transform=fig.transFigure,
        facecolor=ROW_COLORS[disp_label], edgecolor='#888888',
        linewidth=0.8, zorder=2, clip_on=False
    )
    fig.add_artist(box)

    fig.text(
        (box_left + box_right) / 2, (box_bot + box_top) / 2,
        disp_label,
        ha='center', va='center',
        fontsize=20, fontweight='bold', fontfamily='Arial',
        rotation=90, color='#111111', zorder=3
    )

# ── Footnote ───────────────────────────────────────────────────────────────────
fig.text(0.50, 0.018,
         'AW-M-3 replicate excluded from Acute Withdrawal (global outlier across all compartments).',
         ha='center', va='bottom', fontsize=10, color='#666666',
         style='italic', fontfamily='Arial')

pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=180, bbox_inches='tight')
pdf.close()
plt.close(fig)
print(f'Saved: {OUT}')
