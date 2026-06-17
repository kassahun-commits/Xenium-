"""
Membrane AW combined figure  (mirrors make_chromatin_combined.py)
=================================================================
  Top:    Two volcano plots — Increased (red) | Decreased (blue)
          Jitter ±0.20 applied for display only
  Bottom: Venn increased (left) + Venn decreased (right)
          across AW / Intoxication / PA

AW-M-3 excluded (Memb_AW-M-3 not included in AW suffix list).
No Keep/Review filter (Membrane sheet has none).
Output: Membrane_AW_Combined.pdf
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

FILE = '../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Membrane_AW_Combined.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
JITTER      = 0.20

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_NS    = '#EBEBEB'
C_GRAY  = '#AAAAAA'
C_OTHER = '#C8C8C8'

# AW-M-3 excluded
AW_SUFFIXES = ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2']
I_SUFFIXES  = ['I-F-1',  'I-F-2',  'I-F-3',  'I-M-1',  'I-M-2']
PA_SUFFIXES = ['PA-F-1', 'PA-F-2', 'PA-F-3', 'PA-M-1', 'PA-M-2']
N_SUFFIXES  = ['N-F-1',  'N-F-2',  'N-F-3',  'N-M-1',  'N-M-2']
PREFIX      = 'Memb_'

PATCH_IDS  = ['100', '010', '110', '001', '101', '011', '111']
AW_PATCHES = {'100', '110', '101', '111'}

# ── Load data ──────────────────────────────────────────────────────────────────
raw = pd.read_excel(FILE, sheet_name='Membrane')
df  = raw.reset_index(drop=True)

n_cols  = [c for c in df.columns if any(PREFIX + s in str(c) for s in N_SUFFIXES)]
i_cols  = [c for c in df.columns if any(PREFIX + s in str(c) for s in I_SUFFIXES)]
aw_cols = [c for c in df.columns if any(PREFIX + s in str(c) for s in AW_SUFFIXES)]
pa_cols = [c for c in df.columns if any(PREFIX + s in str(c) for s in PA_SUFFIXES)]

print(f'Columns found — N:{len(n_cols)}  I:{len(i_cols)}  AW:{len(aw_cols)}  PA:{len(pa_cols)}')

# ── Stats ──────────────────────────────────────────────────────────────────────
def calc_stats(naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pvs   = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pvs.append(stats.ttest_ind(n, c, equal_var=False)[1]
                   if len(n) >= 2 and len(c) >= 2 else np.nan)
    pvs   = pd.Series(pvs, index=df.index)
    valid = pvs.notna()
    ranks = pvs[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pvs[valid] * valid.sum() / ranks))
    return fc, corr

def get_sig_genes(naive_cols, cond_cols, direction):
    fc, corr = calc_stats(naive_cols, cond_cols)
    mask = fc.notna() & corr.notna() & (corr > CORR_THRESH)
    sig  = mask & (fc > FC_THRESH if direction == 'up' else fc < -FC_THRESH)
    return set(df.loc[sig, 'Gene symbol'].astype(str).tolist())

fc, corrected = calc_stats(n_cols, aw_cols)

aw_up = get_sig_genes(n_cols, aw_cols, 'up')
i_up  = get_sig_genes(n_cols, i_cols,  'up')
pa_up = get_sig_genes(n_cols, pa_cols, 'up')
aw_dn = get_sig_genes(n_cols, aw_cols, 'down')
i_dn  = get_sig_genes(n_cols, i_cols,  'down')
pa_dn = get_sig_genes(n_cols, pa_cols, 'down')

print(f'Membrane AW:  ↑{len(aw_up)}  ↓{len(aw_dn)}  (AW-M-3 excluded)')
print(f'Membrane Intox: ↑{len(i_up)}  ↓{len(i_dn)}')
print(f'Membrane PA:    ↑{len(pa_up)}  ↓{len(pa_dn)}')

# ── Dynamic axis limits ────────────────────────────────────────────────────────
mask    = fc.notna() & corrected.notna()
f_vals  = fc[mask].values
cr_vals = corrected[mask].values

fc_pad  = (f_vals.max() - f_vals.min()) * 0.08
cr_pad  = cr_vals.max() * 0.06
XLIM    = (f_vals.min()  - fc_pad, f_vals.max()  + fc_pad)
YLIM    = (-0.5,                   cr_vals.max()  + cr_pad)

print(f'Axis limits — x:{XLIM}  y:{YLIM}')

# ── Venn helpers ───────────────────────────────────────────────────────────────
def subset_sizes(A, B, C):
    return (len(A-B-C), len(B-A-C), len((A&B)-C),
            len(C-A-B), len((A&C)-B), len((B&C)-A), len(A&B&C))

def make_layout_subsets(A, B, C):
    sA = max(1, len(A)) ** 0.5
    sB = max(1, len(B)) ** 0.5
    sC = max(1, len(C)) ** 0.5
    ov = min(sA, sB, sC) * 0.18
    return (max(0.05, sA - 2*ov), max(0.05, sB - 2*ov), ov,
            max(0.05, sC - 2*ov), ov, ov, ov * 0.4)

def draw_venn(ax, A, B, C, set_labels, highlight_color):
    actual = subset_sizes(A, B, C)
    layout = make_layout_subsets(A, B, C)
    v = venn3(subsets=layout, set_labels=set_labels, ax=ax)

    for pid in PATCH_IDS:
        patch = v.get_patch_by_id(pid)
        if patch:
            if pid in AW_PATCHES:
                patch.set_facecolor(highlight_color); patch.set_alpha(0.75)
            else:
                patch.set_facecolor(C_GRAY); patch.set_alpha(0.30)
            patch.set_edgecolor('none')

    for lbl, val in zip(v.subset_labels, actual):
        if lbl is not None:
            lbl.set_text(str(val)); lbl.set_fontsize(20)
            lbl.set_fontfamily('Arial'); lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    vc = venn3_circles(subsets=layout, ax=ax)
    for circle in (vc.patches if hasattr(vc, 'patches') else vc):
        circle.set_edgecolor('#333333'); circle.set_linewidth(3.5)
        circle.set_facecolor('none')

    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(22); lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold'); lbl.set_color('#111111')

    xl, yl = ax.get_xlim(), ax.get_ylim()
    xc = (xl[0]+xl[1])/2; yc = (yl[0]+yl[1])/2
    xr = (xl[1]-xl[0])*0.70; yr = (yl[1]-yl[0])*0.70
    ax.set_xlim(xc - xr, xc + xr); ax.set_ylim(yc - yr, yc + yr)
    ax.set_axis_off()

# ── Figure ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 16))
fig.patch.set_facecolor('white')

outer = gridspec.GridSpec(2, 1, figure=fig,
                          height_ratios=[1.1, 1.3],
                          hspace=0.18,
                          top=0.93, bottom=0.03,
                          left=0.07, right=0.97)

gs_vp = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], wspace=0.25)
gs_vn = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], wspace=0.45)

ax_vp_up   = fig.add_subplot(gs_vp[0, 0])
ax_vp_dn   = fig.add_subplot(gs_vp[0, 1])
ax_venn_up = fig.add_subplot(gs_vn[0, 0])
ax_venn_dn = fig.add_subplot(gs_vn[0, 1])

# ── Volcano plots ──────────────────────────────────────────────────────────────
np.random.seed(42)
fx   = f_vals + np.random.uniform(-JITTER, JITTER, size=len(f_vals))

up   = (cr_vals > CORR_THRESH) & (f_vals >  FC_THRESH)
down = (cr_vals > CORR_THRESH) & (f_vals < -FC_THRESH)
ns   = ~up & ~down
n_up   = int(up.sum())
n_down = int(down.sum())
n_tot  = int(mask.sum())

def style_vp_ax(ax, hi, lo, c_hi, title, direction, n_hi):
    ax.scatter(fx[ns], cr_vals[ns], s=10, color=C_NS,    alpha=0.40, linewidths=0, rasterized=True)
    ax.scatter(fx[lo], cr_vals[lo], s=18, color=C_OTHER, alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(fx[hi], cr_vals[hi], s=35, color=c_hi,    alpha=0.90, linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlim(XLIM); ax.set_ylim(YLIM)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False);  ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2); ax.spines['bottom'].set_linewidth(1.2)
    ax.set_xlabel('Log$_2$ Fold Change',           fontsize=18, labelpad=8)
    ax.set_ylabel('Corrected p-value\n(−log$_2$)', fontsize=18, labelpad=8)
    ax.set_title(title, fontsize=18, fontweight='bold', pad=10)
    ax.tick_params(labelsize=13)

    dir_sym   = '↑' if direction == 'up' else '↓'
    dir_label = 'Increased' if direction == 'up' else 'Decreased'
    ax.text(0.98, 0.98, f'{dir_sym}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=17, color='#333333', fontweight='bold')
    ax.legend(handles=[
        mpatches.Patch(color=c_hi, label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_NS, label='NS'),
    ], fontsize=12, loc='upper left', framealpha=0.88,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.7)

style_vp_ax(ax_vp_up, up,   down, C_UP,   'Membrane AW vs Naïve — Increased', 'up',   n_up)
style_vp_ax(ax_vp_dn, down, up,   C_DOWN, 'Membrane AW vs Naïve — Decreased', 'down', n_down)

# ── Venn diagrams ──────────────────────────────────────────────────────────────
draw_venn(ax_venn_up, aw_up, i_up, pa_up,
          set_labels=(f'AW Increased\n(n={len(aw_up)})',
                       f'Intox. Increased\n(n={len(i_up)})',
                       f'PA Increased\n(n={len(pa_up)})'),
          highlight_color=C_UP)

draw_venn(ax_venn_dn, aw_dn, i_dn, pa_dn,
          set_labels=(f'AW Decreased\n(n={len(aw_dn)})',
                       f'Intox. Decreased\n(n={len(i_dn)})',
                       f'PA Decreased\n(n={len(pa_dn)})'),
          highlight_color=C_DOWN)

# Panel labels
for ax, letter in zip([ax_vp_up, ax_vp_dn, ax_venn_up, ax_venn_dn], ['A', 'B', 'C', 'D']):
    ax.text(-0.05, 1.05, letter, transform=ax.transAxes,
            fontsize=26, fontweight='bold', fontfamily='Arial',
            va='top', ha='left', color='#111111')

fig.text(0.5, 0.005,
         '|FC| > 0.5, corrected p > 3.3  |  AW-M-3 excluded  |  '
         'Jitter ±0.20 applied for display only; counts are from original data.',
         ha='center', va='bottom', fontsize=8.5,
         color='#666666', style='italic', fontfamily='Arial')

with PdfPages(OUT) as pdf:
    pdf.savefig(fig, dpi=180, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {OUT}')
