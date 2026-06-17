"""
Publication-quality 4×3 volcano plot panel figure.
Rows = compartments (Membrane, Cytosol, Chromatin, Soluble Nuclear)
Cols = conditions (Intoxication, Acute Withdrawal, Protracted Abstinence)
- AW-M-3 excluded from AW only (it is only in the AW group)
- Each panel shows both Increased (red) and Decreased (blue) dots
- Column headers and colored row labels on top/left
- X axis = Log2 Fold Change

Produces THREE PDFs:
  VP_Panel_Figure_v1_jitter40.pdf  — jitter ±0.40 (more spread, for reference)
  VP_Panel_Figure_v2_jitter20.pdf  — jitter ±0.20 (less spread, for reference)
  VP_Panel_Figure_v3_nojitter.pdf  — NO jitter, raw data (journal-safe, recommended)
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

CORR_THRESH = 3.3
FC_THRESH   = 0.5
# Same axes as the Chromatin AW individual VP — applied to all panels
XLIM = (-9.00, 9.74)
YLIM = (-0.3,  27.67)

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_NS    = '#EBEBEB'

ROW_COLORS = {
    'Membrane':        '#E8E4F0',
    'Cytosol':         '#D6EDDC',
    'Soluble Nuclear': '#D4E8F5',
    'Chromatin':       '#F5F0D4',
}

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'Memb_N-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'Cyto_N-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'Chrom_N-', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'Nuc_N-',   'all'),
]

COND_SUFFIXES = {
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2'],   # AW-M-3 excluded
    'Protracted Abstinence': ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
CONDITIONS = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

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

# ── Pre-load all data ──────────────────────────────────────────────────────────
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
        all_data[(ri, ci)] = (fc, corrected, disp_label, len(df))
    print(f'Loaded {disp_label}')

# ── Build one figure per version ───────────────────────────────────────────────
# JITTER=0 means no jitter (raw data, journal-safe)
VERSIONS = [
    ('VP_Panel_Figure_v1_jitter40.pdf', 0.40),
    ('VP_Panel_Figure_v2_jitter20.pdf', 0.20),
    ('VP_Panel_Figure_v3_nojitter.pdf', 0.00),
]

N_ROWS = len(SHEETS)
N_COLS = len(CONDITIONS)
COL_TITLES = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

for OUT, JITTER in VERSIONS:
    fig = plt.figure(figsize=(26, 26))
    fig.patch.set_facecolor('white')

    # left=0.24: axes spine starts at 24% of figure width.
    # y-axis tick labels + rotated title extend ~5-6% left of spine.
    # Colored row-label boxes sit at ~10-16% → clear gap on both sides.
    gs = gridspec.GridSpec(
        N_ROWS, N_COLS,
        figure=fig,
        hspace=0.38,
        wspace=0.32,
        top=0.91, bottom=0.06,
        left=0.24, right=0.97
    )

    axes_grid = {}
    for ri in range(N_ROWS):
        for ci in range(N_COLS):
            fc, corrected, disp_label, n_tot = all_data[(ri, ci)]
            ax = fig.add_subplot(gs[ri, ci])
            axes_grid[(ri, ci)] = ax

            mask = fc.notna() & corrected.notna()
            f    = fc[mask].values
            cr   = corrected[mask].values
            up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
            down = (cr > CORR_THRESH) & (f < -FC_THRESH)
            ns   = ~up & ~down
            n_up   = int(up.sum())
            n_down = int(down.sum())

            # Apply jitter to x positions for display only (0 = no jitter)
            if JITTER > 0:
                np.random.seed(42)
                fx = f + np.random.uniform(-JITTER, JITTER, size=len(f))
            else:
                fx = f   # raw data, no modification

            # Scatter layers
            ax.scatter(fx[ns],   cr[ns],   s=6,  color=C_NS,   alpha=0.35, linewidths=0, rasterized=True)
            ax.scatter(fx[down], cr[down], s=30, color=C_DOWN, alpha=0.85, linewidths=0, rasterized=True)
            ax.scatter(fx[up],   cr[up],   s=30, color=C_UP,   alpha=0.85, linewidths=0, rasterized=True)

            # Threshold lines
            ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.4, alpha=0.7)
            ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.4, alpha=0.7)
            ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.4, alpha=0.7)
            ax.axhline(0,           color='#BBBBBB', linestyle='-',  linewidth=0.7)
            ax.axvline(0,           color='#BBBBBB', linestyle='-',  linewidth=0.7)

            ax.set_xlim(XLIM); ax.set_ylim(YLIM)
            ax.set_facecolor('white')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.0)
            ax.spines['bottom'].set_linewidth(1.0)
            ax.tick_params(labelsize=12)

            ax.set_xlabel('Log$_2$ Fold Change',        fontsize=14, labelpad=5)
            ax.set_ylabel('-Log$_2$ Corrected\np-value', fontsize=14, labelpad=5)

            ax.legend(handles=[
                mpatches.Patch(color=C_UP,   label=f'Increased  (n={n_up})'),
                mpatches.Patch(color=C_DOWN, label=f'Decreased  (n={n_down})'),
                mpatches.Patch(color=C_NS,   label='NS'),
            ], fontsize=10, loc='upper left', framealpha=0.90,
               edgecolor='#CCCCCC', frameon=True, borderpad=0.6,
               handlelength=1.0, labelspacing=0.35)

    # ── Column headers ─────────────────────────────────────────────────────────
    for ci, title in enumerate(COL_TITLES):
        pos = axes_grid[(0, ci)].get_position()
        fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.018, title,
                 ha='center', va='bottom',
                 fontsize=22, fontweight='bold', fontfamily='Arial', color='#111111')

    # ── Row labels — colored strip fully left of y-axis content ───────────────
    for ri, (disp_label, *_) in enumerate(SHEETS):
        pos0 = axes_grid[(ri, 0)].get_position()

        box_right = pos0.x0 - 0.080   # clears tick labels + rotated axis title
        box_left  = box_right - 0.038
        box_bot   = pos0.y0
        box_top   = pos0.y1

        box = FancyBboxPatch(
            (box_left, box_bot),
            box_right - box_left,
            box_top - box_bot,
            boxstyle='round,pad=0.002',
            transform=fig.transFigure,
            facecolor=ROW_COLORS[disp_label],
            edgecolor='#888888',
            linewidth=0.8,
            zorder=2,
            clip_on=False
        )
        fig.add_artist(box)

        fig.text(
            (box_left + box_right) / 2,
            (box_bot  + box_top)  / 2,
            disp_label,
            ha='center', va='center',
            fontsize=20, fontweight='bold', fontfamily='Arial',
            rotation=90, color='#111111',
            zorder=3
        )

    # ── Footnote ───────────────────────────────────────────────────────────────
    fig.text(0.50, 0.022,
             'AW-M-3 replicate excluded from Acute Withdrawal (global outlier across all compartments).',
             ha='center', va='bottom', fontsize=10, color='#666666',
             style='italic', fontfamily='Arial')

    pdf = PdfPages(OUT)
    pdf.savefig(fig, dpi=200, bbox_inches='tight')
    pdf.close()
    plt.close(fig)
    label = f'jitter ±{JITTER}' if JITTER > 0 else 'NO jitter'
    print(f'Saved: {OUT}  ({label})')
