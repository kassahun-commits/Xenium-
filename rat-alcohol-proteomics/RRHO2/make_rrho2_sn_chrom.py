"""
RRHO2: Soluble Nuclear vs Chromatin — one plot per condition (Intox, AW, PA)
================================================================================
Each plot compares the ranked protein lists between the two compartments
for the SAME condition vs Naïve.

  x-axis = Soluble Nuclear  (sign(FC) × corrected_p  vs Naïve)
  y-axis = Chromatin        (sign(FC) × corrected_p  vs Naïve)
  Union of proteins detected in either compartment; missing → score 0

Filters  : Chromatin = Keep+Review  |  Soluble Nuclear = all rows
Outlier  : AW-M-3 excluded from AW
Output   : RRHO2_SN_vs_Chromatin.pdf
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'RRHO2_SN_vs_Chromatin.pdf'

NAIVE_SUFFIXES = ['N-F-1', 'N-F-2', 'N-F-3', 'N-M-1', 'N-M-2']
COND_SUFFIXES  = {
    'Intoxication':          ['I-F-1',  'I-F-2',  'I-F-3',  'I-M-1',  'I-M-2'],
    'Acute Withdrawal':      ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2'],  # AW-M-3 excluded
    'Protracted Abstinence': ['PA-F-1', 'PA-F-2', 'PA-F-3', 'PA-M-1', 'PA-M-2'],
}
CONDITIONS = list(COND_SUFFIXES.keys())

# ── Colormap ───────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'rrho2',
    ['#03007D', '#0033CC', '#0099FF', '#00CCAA',
     '#88DD00', '#FFCC00', '#FF5500', '#880000'],
    N=512
)
CMAP.set_bad('white')

# ── Stats helpers ──────────────────────────────────────────────────────────────
def calc_stats(df, naive_cols, cond_cols):
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
    sc = np.sign(fc) * corr
    sc[fc.isna() | corr.isna()] = np.nan
    return sc

# ── RRHO2 core ─────────────────────────────────────────────────────────────────
def compute_rrho2(score_a, score_b):
    sa = np.where(np.isnan(score_a), 0.0, score_a)
    sb = np.where(np.isnan(score_b), 0.0, score_b)
    n  = len(sa)

    step       = max(1, n // 100)
    thresholds = np.arange(step, n + 1, step)
    nt         = len(thresholds)

    rank_a_up = np.argsort(np.argsort(-sa))
    rank_a_dn = np.argsort(np.argsort( sa))
    rank_b_up = np.argsort(np.argsort(-sb))
    rank_b_dn = np.argsort(np.argsort( sb))

    def one_quadrant(ra, rb):
        a_in    = (ra[np.newaxis, :] < thresholds[:, np.newaxis]).astype(np.int32)
        b_in    = (rb[np.newaxis, :] < thresholds[:, np.newaxis]).astype(np.int32)
        overlap = (a_in @ b_in.T).astype(float)
        K       = thresholds[:, np.newaxis].astype(float)
        n_draw  = thresholds[np.newaxis, :].astype(float)
        pv      = hypergeom.sf(overlap - 1, n, K, n_draw)
        pv      = np.clip(pv, 1e-300, 1.0)
        return -np.log10(pv)

    mat_pp = one_quadrant(rank_a_up, rank_b_up)
    mat_nn = one_quadrant(rank_a_dn, rank_b_dn)
    mat_pn = one_quadrant(rank_a_up, rank_b_dn)
    mat_np = one_quadrant(rank_a_dn, rank_b_up)

    divider     = np.full((nt, 1), np.nan)
    left        = np.concatenate([mat_pp,           divider, np.fliplr(mat_np)],            axis=1)
    right_part  = np.concatenate([np.flipud(mat_pn), divider, np.flipud(np.fliplr(mat_nn))], axis=1)
    divider_row = np.full((1, left.shape[1]), np.nan)
    assembled   = np.concatenate([left, divider_row, right_part], axis=0)
    return assembled, nt

# ── Load Chromatin ─────────────────────────────────────────────────────────────
print('Loading Chromatin (Keep+Review)...')
raw_c  = pd.read_excel(FILE, sheet_name='Chromatin')
df_c   = raw_c[raw_c['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
genes_c = df_c['Gene symbol'].astype(str).values
n_cols_c = [col for col in df_c.columns if 'Chrom_N-' in str(col)]

chrom_scores = {}
for cond, suffixes in COND_SUFFIXES.items():
    cond_cols = [col for col in df_c.columns if any('Chrom_' + s in str(col) for s in suffixes)]
    fc, corr  = calc_stats(df_c, n_cols_c, cond_cols)
    sc        = signed_score(fc, corr)
    s         = pd.Series(sc.values, index=genes_c)
    chrom_scores[cond] = s[~s.index.duplicated(keep='first')]
    print(f'  Chromatin {cond}: {s.notna().sum()} valid scores')

# ── Load Soluble Nuclear ───────────────────────────────────────────────────────
print('Loading Soluble Nuclear (all rows)...')
raw_n   = pd.read_excel(FILE, sheet_name='Soluble nuclear')
df_n    = raw_n.copy().reset_index(drop=True)
genes_n = df_n['Gene symbol'].astype(str).values
n_cols_n = [col for col in df_n.columns if 'Nuc_N-' in str(col)]

nuc_scores = {}
for cond, suffixes in COND_SUFFIXES.items():
    cond_cols = [col for col in df_n.columns if any('Nuc_' + s in str(col) for s in suffixes)]
    fc, corr  = calc_stats(df_n, n_cols_n, cond_cols)
    sc        = signed_score(fc, corr)
    s         = pd.Series(sc.values, index=genes_n)
    nuc_scores[cond] = s[~s.index.duplicated(keep='first')]
    print(f'  SN {cond}: {s.notna().sum()} valid scores')

# ── Run RRHO2 for each condition ───────────────────────────────────────────────
print('\nRunning RRHO2...')
results     = {}
vmax_global = 0.0

for cond in CONDITIONS:
    sn_s   = nuc_scores[cond]
    ch_s   = chrom_scores[cond]

    all_genes = sorted(set(sn_s.index) | set(ch_s.index))
    n         = len(all_genes)

    sa = np.array([sn_s.get(g, np.nan)  for g in all_genes], dtype=float)
    sb = np.array([ch_s.get(g, np.nan)  for g in all_genes], dtype=float)

    mat, nt     = compute_rrho2(sa, sb)
    vmax_global = max(vmax_global, np.nanmax(mat))
    results[cond] = (mat, nt, n)
    print(f'  {cond}: n={n}, nt={nt}, max={np.nanmax(mat):.1f}')

print(f'\nGlobal max = {vmax_global:.1f}')

# ── Figure: 1 row × 3 columns ─────────────────────────────────────────────────
N_COLS = len(CONDITIONS)

FIG_W = 4.2 * N_COLS + 1.8    # extra for colorbar
FIG_H = 5.2                    # enough for labels + quadrant key

fig = plt.figure(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')

gs = GridSpec(1, N_COLS,
              figure=fig,
              left=0.06, right=0.88,
              top=0.88,  bottom=0.30,
              hspace=0.0, wspace=0.45)

axes = [fig.add_subplot(gs[0, ci]) for ci in range(N_COLS)]

im_ref = None
for ci, cond in enumerate(CONDITIONS):
    ax         = axes[ci]
    mat, nt, n = results[cond]

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
        spine.set_linewidth(0.7)
        spine.set_color('#888888')

    ax.set_title(cond, fontsize=12, fontweight='bold', pad=8)
    ax.set_xlabel('Soluble Nuclear\n(fold change vs Naïve)',
                  fontsize=9.5, labelpad=6, color='#222222')
    ax.set_ylabel('Chromatin\n(fold change vs Naïve)',
                  fontsize=9.5, labelpad=6, color='#222222')

# ── Colorbar ───────────────────────────────────────────────────────────────────
cbar_ax = fig.add_axes([0.905, 0.30, 0.020, 0.58])
cb = fig.colorbar(im_ref, cax=cbar_ax)
cb.set_label('−log₁₀(p-value)', fontsize=10, labelpad=8)
cb.ax.tick_params(labelsize=9)

# ── Quadrant key ───────────────────────────────────────────────────────────────
key_ax = fig.add_axes([0.28, 0.01, 0.26, 0.22])
key_ax.set_facecolor('white')

hot_color  = '#CC2200'
cold_color = '#0044BB'

boxes = [
    (0.05, 0.52, 0.42, 0.44, hot_color,  'Both ↑\n(SN & Chromatin)',  'white'),
    (0.53, 0.52, 0.42, 0.44, cold_color, 'SN ↓  Chrom ↑\n(discordant)', 'white'),
    (0.05, 0.04, 0.42, 0.44, cold_color, 'SN ↑  Chrom ↓\n(discordant)', 'white'),
    (0.53, 0.04, 0.42, 0.44, hot_color,  'Both ↓\n(SN & Chromatin)',  'white'),
]

for x0, y0, w, h, color, label, fc in boxes:
    patch = FancyBboxPatch((x0, y0), w, h,
                           boxstyle='round,pad=0.02',
                           facecolor=color, edgecolor='white', linewidth=1.5,
                           transform=key_ax.transAxes, clip_on=False)
    key_ax.add_patch(patch)
    key_ax.text(x0 + w/2, y0 + h/2, label,
                transform=key_ax.transAxes,
                fontsize=7.5, color=fc, ha='center', va='center',
                fontweight='bold', multialignment='center')

# Axis arrows
key_ax.annotate('', xy=(0.97, 0.50), xytext=(0.03, 0.50),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))
key_ax.annotate('', xy=(0.50, 0.03), xytext=(0.50, 0.97),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5))

key_ax.text(0.50, -0.10, 'Soluble Nuclear  (↑ left  →  ↓ right)',
            transform=key_ax.transAxes, fontsize=8, ha='center', va='top', color='#333333')
key_ax.text(-0.08, 0.50, 'Chromatin\n(↑ top → ↓ bottom)',
            transform=key_ax.transAxes, fontsize=8, ha='right', va='center',
            color='#333333', multialignment='center')

key_ax.set_xlim(0, 1); key_ax.set_ylim(0, 1); key_ax.axis('off')
key_ax.text(0.50, 1.10, 'Quadrant Guide',
            transform=key_ax.transAxes, fontsize=10, ha='center', va='bottom',
            fontweight='bold', color='#222222')

with PdfPages(OUT) as pdf:
    pdf.savefig(fig, dpi=180, bbox_inches='tight')
plt.close()

print(f'\nSaved: {OUT}')
