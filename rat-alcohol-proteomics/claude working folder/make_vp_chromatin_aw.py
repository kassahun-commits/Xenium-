"""
Chromatin AW volcano plot with jitter on x-axis to prevent dot stacking artifact.
"""
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'VP_Chromatin_AW_jitter.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
JITTER      = 0.04   # small random nudge on x-axis (in FC units)

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_OTHER = '#C8C8C8'
C_NS    = '#EBEBEB'

np.random.seed(42)

# ── Load data ──────────────────────────────────────────────────────────────────
df = pd.read_excel(FILE, sheet_name='Chromatin', header=0)
df = df[df['Filter'].isin(['Keep','Review'])].reset_index(drop=True)

naive_cols = [c for c in df.columns if 'Chrom_N-' in str(c)]
aw_cols    = [c for c in df.columns
              if any('Chrom_'+s in str(c)
                     for s in ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'])]

naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
cond  = df[aw_cols].apply(pd.to_numeric, errors='coerce')
fc    = cond.mean(axis=1) - naive.mean(axis=1)

pv = []
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

mask = fc.notna() & corr.notna()
f, cr = fc[mask].values, corr[mask].values

up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
down = (cr > CORR_THRESH) & (f < -FC_THRESH)
ns   = ~up & ~down

n_up   = up.sum()
n_down = down.sum()
n_tot  = mask.sum()

# Add jitter only to overlapping x positions
jitter_x = f + np.random.uniform(-JITTER, JITTER, size=len(f))

# ── Plot ───────────────────────────────────────────────────────────────────────
def make_panel(ax, highlight, title):
    hi, lo, c_hi = (up, down, C_UP) if highlight == 'up' else (down, up, C_DOWN)

    ax.scatter(jitter_x[ns],  cr[ns],  s=12, color=C_NS,    alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(jitter_x[lo],  cr[lo],  s=22, color=C_OTHER, alpha=0.50, linewidths=0, rasterized=True)
    ax.scatter(jitter_x[hi],  cr[hi],  s=40, color=c_hi,    alpha=0.90, linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlabel('Fold Change',                    fontsize=20, labelpad=8)
    ax.set_ylabel('Corrected p-value\n(−log₂)',     fontsize=20, labelpad=8)
    ax.set_title(title, fontsize=21, fontweight='bold', pad=12)
    ax.tick_params(labelsize=15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_hi = int(hi.sum())
    direction = '↑' if highlight == 'up' else '↓'
    ax.text(0.98, 0.98, f'{direction}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=18, color='#333333', fontweight='bold')

    dir_label = 'Increased' if highlight == 'up' else 'Decreased'
    ax.legend(handles=[
        mpatches.Patch(color=c_hi,    label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_OTHER, label='Other sig.'),
        mpatches.Patch(color=C_NS,    label='NS'),
    ], fontsize=14, loc='upper left', framealpha=0.88,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.8)

with PdfPages(OUT) as pdf:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.patch.set_facecolor('white')

    make_panel(axes[0], 'up',
               'Chromatin — AW vs Naïve\nIncreased')
    make_panel(axes[1], 'down',
               'Chromatin — AW vs Naïve\nDecreased')

    fig.text(0.5, 0.01,
             f'Thresholds: |FC| > {FC_THRESH}  |  Corrected p > {CORR_THRESH}  '
             f'|  Jitter = ±{JITTER} FC units (x-axis only)',
             ha='center', fontsize=11, color='#666666', style='italic')

    plt.tight_layout(pad=2.0, rect=[0, 0.04, 1, 1])
    pdf.savefig(fig, dpi=180, bbox_inches='tight')
    plt.close(fig)

print(f'Saved: {OUT}')
print(f'↑{n_up}  ↓{n_down}  total={n_tot}')
