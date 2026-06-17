"""
Publication-quality 4×3 volcano plot panel figure.
Rows = compartments (Membrane, Cytosol, Chromatin, Soluble Nuclear)
Cols = conditions (Intoxication, Acute Withdrawal, Protracted Abstinence)
- AW-M-3 excluded from AW
- Each panel shows both Increased (red) and Decreased (blue)
- Column headers and colored row labels on top/left

Produces FOUR PDFs:
  VP_Panel_Figure_v1_jitter40.pdf  — jitter ±0.40 (more spread)
  VP_Panel_Figure_v2_jitter20.pdf  — jitter ±0.20 (less spread)
  VP_Panel_Figure_v3_nojitter.pdf  — NO jitter, raw data (journal-safe)
  VP_Panel_Figure_v4_bigfont.pdf   — no jitter, large fonts (~4× labels, ~3× VP text)
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
XLIM = (-9.00, 9.74)
YLIM = (-0.3,  27.67)

C_UP  = '#E8305A'
C_DOWN = '#2B7FD4'
C_NS   = '#EBEBEB'

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
COL_TITLES  = CONDITIONS

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
        all_data[(ri, ci)] = (fc, corrected, disp_label)
    print(f'Loaded {disp_label}')

# ── Version configs ────────────────────────────────────────────────────────────
# Each dict controls one output PDF.
# fs_hdr   = compartment + group label font size (the big outer labels)
# fs_row   = row-label box text font size
# fs_ax    = axis title font size inside each VP
# fs_tick  = axis tick number font size inside each VP
# fs_leg   = legend font size inside each VP
# gs_left  = gridspec left margin (must leave room for y-axis text + row boxes)
# gs_top   = gridspec top margin (must leave room for column headers)
# gs_bot   = gridspec bottom margin (must leave room for x-axis labels)
# hspace   = vertical gap between rows
# wspace   = horizontal gap between columns
# box_w    = width of colored row-label strip (figure fraction)
# box_gap  = clearance between row-label strip right edge and axes spine
# dot_ns   = scatter dot size for NS proteins
# dot_sig  = scatter dot size for significant proteins

VERSIONS = [
    dict(
        out    = 'VP_Panel_Figure_v1_jitter40.pdf',
        jitter = 0.40,
        fs_hdr = 22,  fs_row = 20,
        fs_ax  = 14,  fs_tick= 12,  fs_leg = 10,
        gs_left= 0.24, gs_top= 0.91, gs_bot= 0.06,
        hspace = 0.38, wspace= 0.32,
        box_w  = 0.038, box_gap= 0.080,
        dot_ns = 6,   dot_sig= 30,
    ),
    dict(
        out    = 'VP_Panel_Figure_v2_jitter20.pdf',
        jitter = 0.20,
        fs_hdr = 22,  fs_row = 20,
        fs_ax  = 14,  fs_tick= 12,  fs_leg = 10,
        gs_left= 0.24, gs_top= 0.91, gs_bot= 0.06,
        hspace = 0.38, wspace= 0.32,
        box_w  = 0.038, box_gap= 0.080,
        dot_ns = 6,   dot_sig= 30,
    ),
    dict(
        out    = 'VP_Panel_Figure_v3_nojitter.pdf',
        jitter = 0.00,
        fs_hdr = 22,  fs_row = 20,
        fs_ax  = 14,  fs_tick= 12,  fs_leg = 10,
        gs_left= 0.24, gs_top= 0.91, gs_bot= 0.06,
        hspace = 0.38, wspace= 0.32,
        box_w  = 0.038, box_gap= 0.080,
        dot_ns = 6,   dot_sig= 30,
    ),
    dict(
        out    = 'VP_Panel_Figure_v4_bigfont.pdf',
        jitter = 0.00,
        # ~4× outer labels, ~3× VP text
        fs_hdr = 88,  fs_row = 80,
        fs_ax  = 42,  fs_tick= 36,  fs_leg = 30,
        # Wider margins to fit larger axis text and row-label boxes
        gs_left= 0.38, gs_top= 0.88, gs_bot= 0.13,
        hspace = 0.60, wspace= 0.50,
        box_w  = 0.065, box_gap= 0.110,
        dot_ns = 18,  dot_sig= 80,
    ),
]

N_ROWS = len(SHEETS)
N_COLS = len(CONDITIONS)

# ── Build each version ─────────────────────────────────────────────────────────
for cfg in VERSIONS:
    fig = plt.figure(figsize=(26, 26))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(
        N_ROWS, N_COLS,
        figure  = fig,
        hspace  = cfg['hspace'],
        wspace  = cfg['wspace'],
        top     = cfg['gs_top'],
        bottom  = cfg['gs_bot'],
        left    = cfg['gs_left'],
        right   = 0.97,
    )

    # ── Volcano panels ─────────────────────────────────────────────────────────
    axes_grid = {}
    for ri in range(N_ROWS):
        for ci in range(N_COLS):
            fc, corrected, disp_label = all_data[(ri, ci)]
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

            # x-axis jitter for display only (0 = raw data, no change)
            if cfg['jitter'] > 0:
                np.random.seed(42)
                fx = f + np.random.uniform(-cfg['jitter'], cfg['jitter'], size=len(f))
            else:
                fx = f

            ax.scatter(fx[ns],   cr[ns],   s=cfg['dot_ns'],  color=C_NS,
                       alpha=0.35, linewidths=0, rasterized=True)
            ax.scatter(fx[down], cr[down], s=cfg['dot_sig'], color=C_DOWN,
                       alpha=0.85, linewidths=0, rasterized=True)
            ax.scatter(fx[up],   cr[up],   s=cfg['dot_sig'], color=C_UP,
                       alpha=0.85, linewidths=0, rasterized=True)

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
            ax.tick_params(labelsize=cfg['fs_tick'])

            ax.set_xlabel('Log$_2$ Fold Change',         fontsize=cfg['fs_ax'], labelpad=5)
            ax.set_ylabel('-Log$_2$ Corrected\np-value', fontsize=cfg['fs_ax'], labelpad=5)

            ax.legend(handles=[
                mpatches.Patch(color=C_UP,   label=f'Increased  (n={n_up})'),
                mpatches.Patch(color=C_DOWN, label=f'Decreased  (n={n_down})'),
                mpatches.Patch(color=C_NS,   label='NS'),
            ], fontsize=cfg['fs_leg'], loc='upper left', framealpha=0.90,
               edgecolor='#CCCCCC', frameon=True, borderpad=0.6,
               handlelength=1.0, labelspacing=0.35)

    # ── Column headers ─────────────────────────────────────────────────────────
    for ci, title in enumerate(COL_TITLES):
        pos = axes_grid[(0, ci)].get_position()
        fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.018, title,
                 ha='center', va='bottom',
                 fontsize=cfg['fs_hdr'], fontweight='bold',
                 fontfamily='Arial', color='#111111')

    # ── Row labels — colored strip fully left of y-axis content ───────────────
    for ri, (disp_label, *_) in enumerate(SHEETS):
        pos0 = axes_grid[(ri, 0)].get_position()

        box_right = pos0.x0 - cfg['box_gap']
        box_left  = box_right - cfg['box_w']
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
            zorder=2, clip_on=False
        )
        fig.add_artist(box)

        fig.text(
            (box_left + box_right) / 2,
            (box_bot  + box_top)  / 2,
            disp_label,
            ha='center', va='center',
            fontsize=cfg['fs_row'], fontweight='bold',
            fontfamily='Arial', rotation=90, color='#111111',
            zorder=3
        )

    # ── Footnote ───────────────────────────────────────────────────────────────
    fig.text(0.50, 0.018,
             'AW-M-3 replicate excluded from Acute Withdrawal (global outlier across all compartments).',
             ha='center', va='bottom', fontsize=max(10, cfg['fs_tick']//2),
             color='#666666', style='italic', fontfamily='Arial')

    pdf = PdfPages(cfg['out'])
    pdf.savefig(fig, dpi=200, bbox_inches='tight')
    pdf.close()
    plt.close(fig)
    jlabel = f"jitter ±{cfg['jitter']}" if cfg['jitter'] > 0 else 'NO jitter'
    print(f"Saved: {cfg['out']}  ({jlabel}, hdr={cfg['fs_hdr']}pt, ax={cfg['fs_ax']}pt)")
