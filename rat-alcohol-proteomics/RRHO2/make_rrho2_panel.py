"""
RRHO2 Panel  –  4 compartments × 3 pairwise condition comparisons (12 plots)
================================================================================
Layout  : 3 rows (pairwise comparisons)  ×  4 cols (compartments)
           Membrane | Cytosol | Chromatin | Soluble Nuclear
           Row 1 : Intox  vs AW
           Row 2 : Intox  vs PA
           Row 3 : AW     vs PA

Ranking  : sign(FC) × corrected_p   (signed metric from Welch's t-test vs Naïve)
           corrected_p = −log2(min(1, p × N / rank))   [higher = more significant]
           Missing score (not detected) → treated as 0 (neutral)

Union    : proteins detected in EITHER condition in that compartment
Outlier  : AW-M-3 excluded from all AW comparisons
Filter   : Chromatin → Keep+Review;  all others → all rows

Output   : RRHO2_Panel.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'RRHO2_Panel.pdf'

# ── Compartment config ─────────────────────────────────────────────────────────
#   (display_label, sheet_name, col_prefix, naive_pattern, filter_mode)
SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'Memb_N-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'Cyto_N-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'Chrom_N-', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'Nuc_N-',   'all'),
]

NAIVE_SUFFIXES = ['N-F-1', 'N-F-2', 'N-F-3', 'N-M-1', 'N-M-2']

COND_SUFFIXES = {
    'Intox': ['I-F-1',  'I-F-2',  'I-F-3',  'I-M-1',  'I-M-2'],
    'AW':    ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2'],   # AW-M-3 excluded
    'PA':    ['PA-F-1', 'PA-F-2', 'PA-F-3', 'PA-M-1', 'PA-M-2'],
}

COND_LABELS = {
    'Intox': 'Intoxication',
    'AW':    'Acute Withdrawal',
    'PA':    'Prot. Abstinence',
}

# 3 pairwise comparisons  (x-condition, y-condition)
COMPARISONS = [
    ('Intox', 'AW'),
    ('Intox', 'PA'),
    ('AW',    'PA'),
]

# ── Colormap (dark-blue → blue → cyan → green → yellow → orange → red) ────────
CMAP = LinearSegmentedColormap.from_list(
    'rrho2',
    ['#03007D', '#0033CC', '#0099FF', '#00CCAA',
     '#88DD00', '#FFCC00', '#FF5500', '#880000'],
    N=512
)
CMAP.set_bad('white')   # NaN → white (used for the divider gap)

# ── Statistics ─────────────────────────────────────────────────────────────────
def calc_stats(df, naive_cols, cond_cols):
    """Returns (fc, corrected_p) as pd.Series aligned to df.index."""
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pvs   = []
    for i in range(len(df)):
        nv = naive.iloc[i].dropna().values
        cv = cond.iloc[i].dropna().values
        pvs.append(stats.ttest_ind(nv, cv, equal_var=False)[1]
                   if len(nv) >= 2 and len(cv) >= 2 else np.nan)
    pvs   = pd.Series(pvs, index=df.index)
    valid = pvs.notna()
    ranks = pvs[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1.0, pvs[valid] * valid.sum() / ranks))
    return fc, corr

def signed_score(fc, corr):
    """sign(FC) × corrected_p;  NaN if either input is NaN."""
    sc = np.sign(fc) * corr
    sc[fc.isna() | corr.isna()] = np.nan
    return sc

# ── RRHO2 core ─────────────────────────────────────────────────────────────────
def compute_rrho2(score_a, score_b):
    """
    Parameters
    ----------
    score_a, score_b : 1-D arrays (same length)
        Signed ranking scores.  NaN treated as 0 (neutral).

    Returns
    -------
    assembled : ndarray  shape (2*nt+1, 2*nt+1)   (+1 for white divider row/col)
        -log10 hypergeometric p-values arranged in 4 quadrants:
          top-left  = both upregulated   (concordant ↑)
          bot-right = both downregulated (concordant ↓)
          top-right = down A / up B      (discordant)
          bot-left  = up A / down B      (discordant)
    nt : int  (number of threshold steps per quadrant)
    """
    sa = np.where(np.isnan(score_a), 0.0, score_a)
    sb = np.where(np.isnan(score_b), 0.0, score_b)
    n  = len(sa)

    step       = max(1, n // 100)          # ~100 steps → smooth heatmap
    thresholds = np.arange(step, n + 1, step)
    nt         = len(thresholds)

    # Rank each gene in descending (up) and ascending (down) order
    rank_a_up = np.argsort(np.argsort(-sa))   # 0 = most positive
    rank_a_dn = np.argsort(np.argsort( sa))   # 0 = most negative
    rank_b_up = np.argsort(np.argsort(-sb))
    rank_b_dn = np.argsort(np.argsort( sb))

    def one_quadrant(ra, rb):
        # a_in[i,g] = 1 if gene g is in top thresholds[i] by ranking ra
        a_in = (ra[np.newaxis, :] < thresholds[:, np.newaxis]).astype(np.int32)  # (nt, n)
        b_in = (rb[np.newaxis, :] < thresholds[:, np.newaxis]).astype(np.int32)  # (nt, n)
        overlap = (a_in @ b_in.T).astype(float)   # (nt, nt)  — fast BLAS multiply

        K      = thresholds[:, np.newaxis].astype(float)   # set size in A
        n_draw = thresholds[np.newaxis, :].astype(float)   # set size in B

        pv = hypergeom.sf(overlap - 1, n, K, n_draw)
        pv = np.clip(pv, 1e-300, 1.0)
        return -np.log10(pv)

    mat_pp = one_quadrant(rank_a_up, rank_b_up)   # both up
    mat_nn = one_quadrant(rank_a_dn, rank_b_dn)   # both down
    mat_pn = one_quadrant(rank_a_up, rank_b_dn)   # up A, down B
    mat_np = one_quadrant(rank_a_dn, rank_b_up)   # down A, up B

    # ── Assemble ──────────────────────────────────────────────────────────────
    # x-axis: left = most upregulated in A  →  right = most downregulated in A
    # y-axis: top  = most upregulated in B  →  bottom = most downregulated in B
    #
    # mat_pp[i,j]: i=0 = smallest up-threshold for A, j=0 for B
    #   → most extreme corner should be at top-left  → place directly
    # mat_np[i,j]: i=0 = most-down-A, j=0 = most-up-B
    #   → extreme corner at top-right  → flip columns
    # mat_pn[i,j]: i=0 = most-up-A, j=0 = most-down-B
    #   → extreme corner at bottom-left → flip rows
    # mat_nn[i,j]: i=0 = most-down-A, j=0 = most-down-B
    #   → extreme corner at bottom-right → flip both
    #
    # We insert a 1-pixel white divider row and column between the quadrants.

    divider = np.full((nt, 1), np.nan)         # vertical divider column
    left  = np.concatenate([mat_pp,                     divider, np.fliplr(mat_np)],          axis=1)
    right = np.concatenate([np.flipud(mat_pn),          divider, np.flipud(np.fliplr(mat_nn))], axis=1)

    divider_row = np.full((1, left.shape[1]), np.nan)  # horizontal divider row
    assembled = np.concatenate([left, divider_row, right], axis=0)

    return assembled, nt

# ── Load data & compute scores ─────────────────────────────────────────────────
print('Loading data and computing scores...')
comp_scores = {}   # {disp_label: {cond: pd.Series indexed by gene symbol}}

for disp_label, sheet, prefix, n_pat, fmode in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    genes      = df['Gene symbol'].astype(str).values
    naive_cols = [c for c in df.columns if n_pat in str(c)]

    cond_data = {}
    for cond_name, suffixes in COND_SUFFIXES.items():
        cond_cols = [c for c in df.columns if any(prefix + s in str(c) for s in suffixes)]
        fc, corr  = calc_stats(df, naive_cols, cond_cols)
        sc        = signed_score(fc, corr)
        # Build series and deduplicate gene symbols (keep first non-NaN, else first)
        s = pd.Series(sc.values, index=genes)
        s = s[~s.index.duplicated(keep='first')]
        cond_data[cond_name] = s

    comp_scores[disp_label] = cond_data
    n_valid = {k: int(v.notna().sum()) for k, v in cond_data.items()}
    print(f'  {disp_label}: {len(genes)} proteins  |  valid scores: {n_valid}')

# ── Run RRHO2 for all 12 combinations ─────────────────────────────────────────
print('\nRunning RRHO2 (12 pairwise × 4 quadrants)...')
rrho2_results = {}   # {(comp_idx, disp_label): (assembled_matrix, nt, n)}
vmax_global = 0.0

for ci, (cond_a, cond_b) in enumerate(COMPARISONS):
    for disp_label, *_ in SHEETS:
        sa_series = comp_scores[disp_label][cond_a]
        sb_series = comp_scores[disp_label][cond_b]

        # Union gene list; missing score → 0 (neutral)
        all_genes = sorted(set(sa_series.index) | set(sb_series.index))
        n = len(all_genes)

        sa = np.array([sa_series.get(g, np.nan) for g in all_genes], dtype=float)
        sb = np.array([sb_series.get(g, np.nan) for g in all_genes], dtype=float)

        mat, nt = compute_rrho2(sa, sb)
        vmax_global = max(vmax_global, np.nanmax(mat))
        rrho2_results[(ci, disp_label)] = (mat, nt, n)
        print(f'  {disp_label:16s}  {cond_a} vs {cond_b}: '
              f'n={n}, nt={nt}, max={np.nanmax(mat):.1f}')

print(f'\nGlobal max -log10(p) = {vmax_global:.1f}')

# ── Figure layout ──────────────────────────────────────────────────────────────
N_ROWS = len(COMPARISONS)
N_COLS = len(SHEETS)

# Extra left margin for row labels, right margin for colorbar
# Extra bottom margin for quadrant legend
FIG_W = 3.6 * N_COLS + 2.8
FIG_H = 3.6 * N_ROWS + 2.8   # extra height for quadrant key below

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')

# GridSpec: leave room on left for row labels, right for colorbar, bottom for key
from matplotlib.gridspec import GridSpec
gs = GridSpec(N_ROWS, N_COLS,
              figure=fig,
              left=0.10, right=0.88,
              top=0.95,  bottom=0.22,
              hspace=0.55, wspace=0.45)

axes = np.empty((N_ROWS, N_COLS), dtype=object)
for ri in range(N_ROWS):
    for ci in range(N_COLS):
        axes[ri, ci] = fig.add_subplot(gs[ri, ci])

# ── Column headers (compartment names) ────────────────────────────────────────
for ci, (disp_label, *_) in enumerate(SHEETS):
    axes[0, ci].set_title(disp_label, fontsize=13, fontweight='bold', pad=9)

# ── Row labels (comparison descriptions) on left side ─────────────────────────
ROW_LABELS = [
    f'{COND_LABELS[ca]}\nvs\n{COND_LABELS[cb]}'
    for ca, cb in COMPARISONS
]
# Compute y-positions from axes positions (set after drawing)
# We'll add them after imshow so the axes are positioned

im_ref = None

for ri, (cond_a, cond_b) in enumerate(COMPARISONS):
    for ci, (disp_label, *_) in enumerate(SHEETS):
        ax = axes[ri, ci]
        mat, nt, n = rrho2_results[(ri, disp_label)]

        im = ax.imshow(mat,
                       aspect='auto',
                       cmap=CMAP,
                       vmin=0, vmax=vmax_global,
                       origin='upper',
                       interpolation='bilinear')
        if im_ref is None:
            im_ref = im

        ax.set_xticks([])
        ax.set_yticks([])

        for spine in ax.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color('#888888')

        # X and Y axis labels on EVERY subplot so each plot is self-contained
        ax.set_xlabel(COND_LABELS[cond_a], fontsize=8, labelpad=4, color='#222222')
        ax.set_ylabel(COND_LABELS[cond_b], fontsize=8, labelpad=4, color='#222222')

# (row labels handled per-subplot below via set_xlabel / set_ylabel)

# ── Colorbar ───────────────────────────────────────────────────────────────────
cbar_ax = fig.add_axes([0.905, 0.22, 0.018, 0.72])
cb = fig.colorbar(im_ref, cax=cbar_ax)
cb.set_label('−log₁₀(p-value)', fontsize=10, labelpad=8)
cb.ax.tick_params(labelsize=9)

# ── Quadrant key diagram ────────────────────────────────────────────────────────
# Placed below the main grid, centered
# Shows a 2×2 schematic of what each quadrant represents
key_ax = fig.add_axes([0.30, 0.01, 0.22, 0.18])   # [left, bottom, width, height]
key_ax.set_facecolor('white')

# Draw 4 colored boxes
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

hot_color  = '#CC2200'   # concordant (both up or both down)
cold_color = '#0044BB'   # discordant
mid_color  = '#DDDDDD'   # border

boxes = [
    # (x0, y0, w, h, color, label, fontcolor)
    (0.05, 0.55, 0.40, 0.40, hot_color,  'Both ↑\nin X & Y',  'white'),   # top-left
    (0.55, 0.55, 0.40, 0.40, cold_color, 'X ↓  Y ↑\n(discordant)', 'white'),  # top-right
    (0.05, 0.05, 0.40, 0.40, cold_color, 'X ↑  Y ↓\n(discordant)', 'white'),  # bottom-left
    (0.55, 0.05, 0.40, 0.40, hot_color,  'Both ↓\nin X & Y',  'white'),   # bottom-right
]

for x0, y0, w, h, color, label, fc in boxes:
    patch = FancyBboxPatch((x0, y0), w, h,
                           boxstyle='round,pad=0.02',
                           facecolor=color, edgecolor='white', linewidth=1.5,
                           transform=key_ax.transAxes, clip_on=False)
    key_ax.add_patch(patch)
    key_ax.text(x0 + w/2, y0 + h/2, label,
                transform=key_ax.transAxes,
                fontsize=8.5, color=fc, ha='center', va='center',
                fontweight='bold', multialignment='center')

# Axis arrows and labels
key_ax.annotate('', xy=(0.98, 0.50), xytext=(0.02, 0.50),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))
key_ax.annotate('', xy=(0.50, 0.02), xytext=(0.50, 0.98),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))

key_ax.text(0.50, -0.08, 'X-axis condition  (↑ left  →  ↓ right)',
            transform=key_ax.transAxes, fontsize=8, ha='center', va='top',
            color='#333333')
key_ax.text(-0.08, 0.50, 'Y-axis\ncondition\n(↑ top\n↓ bottom)',
            transform=key_ax.transAxes, fontsize=8, ha='right', va='center',
            color='#333333', multialignment='center')

key_ax.set_xlim(0, 1)
key_ax.set_ylim(0, 1)
key_ax.axis('off')

# Title for key
key_ax.text(0.50, 1.08, 'Quadrant Guide',
            transform=key_ax.transAxes, fontsize=10, ha='center', va='bottom',
            fontweight='bold', color='#222222')

with PdfPages(OUT) as pdf:
    pdf.savefig(fig, dpi=180, bbox_inches='tight')
plt.close()

print(f'\nSaved: {OUT}')
