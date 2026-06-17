"""
All 12 volcano plots (4 compartments × 3 conditions: Intox, AW, PA vs Naïve)
AW-M-3 excluded from AW comparisons.
Fixed axes: x = ±7 (FC), y = 0–20 (corrected p) for all plots.
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
OUT  = 'VP_All12_noAWM3_FixedAxes.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5
XLIM = (-7, 7)
YLIM = (0, 20)

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_OTHER = '#C8C8C8'
C_NS    = '#EBEBEB'

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'Memb_N-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'Cyto_N-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'Chrom_N-', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'Nuc_N-',   'all'),
]

# AW-M-3 excluded from AW; all 5 replicates used for Intox and PA
COND_SUFFIXES = {
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2'],     # AW-M-3 excluded
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

def draw_volcano(ax, fc, corrected, highlight, title):
    mask = fc.notna() & corrected.notna()
    f    = np.clip(fc[mask].values,        XLIM[0], XLIM[1])
    cr   = np.clip(corrected[mask].values, YLIM[0], YLIM[1])
    fc_r = fc[mask].values
    cr_r = corrected[mask].values

    up   = (cr_r > CORR_THRESH) & (fc_r >  FC_THRESH)
    down = (cr_r > CORR_THRESH) & (fc_r < -FC_THRESH)
    ns   = ~up & ~down
    hi, lo, c_hi = (up, down, C_UP) if highlight == 'up' else (down, up, C_DOWN)

    ax.scatter(f[ns], cr[ns], s=10, color=C_NS,    alpha=0.40, linewidths=0, rasterized=True)
    ax.scatter(f[lo], cr[lo], s=18, color=C_OTHER, alpha=0.50, linewidths=0, rasterized=True)
    ax.scatter(f[hi], cr[hi], s=35, color=c_hi,    alpha=0.90, linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlim(XLIM); ax.set_ylim(YLIM)
    ax.set_xlabel('Fold Change',                fontsize=18, labelpad=8)
    ax.set_ylabel('Corrected p-value\n(−log₂)', fontsize=18, labelpad=8)
    ax.set_title(title, fontsize=19, fontweight='bold', pad=10)
    ax.tick_params(labelsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_hi  = int(hi.sum())
    n_tot = int(mask.sum())
    direction = '↑' if highlight == 'up' else '↓'
    ax.text(0.98, 0.98, f'{direction}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=17, color='#333333', fontweight='bold')

    dir_label = 'Increased' if highlight == 'up' else 'Decreased'
    ax.legend(handles=[
        mpatches.Patch(color=c_hi,    label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_OTHER, label='Other sig.'),
        mpatches.Patch(color=C_NS,    label='NS'),
    ], fontsize=12, loc='upper left', framealpha=0.88,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.8)

with PdfPages(OUT) as pdf:
    for disp_label, sheet, prefix, n_pat, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=sheet)
        if fmode == 'keep_review' and 'Filter' in raw.columns:
            df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        naive_cols = [c for c in df.columns if n_pat in str(c)]

        for cond_name in CONDITIONS:
            suffixes   = COND_SUFFIXES[cond_name]
            cond_cols  = [c for c in df.columns if any(prefix+s in str(c) for s in suffixes)]
            fc, corrected = calc_stats(df, naive_cols, cond_cols)

            aw_note = '  (AW-M-3 excl.)' if cond_name == 'Acute Withdrawal' else ''

            fig, axes = plt.subplots(1, 2, figsize=(16, 6))
            fig.patch.set_facecolor('white')

            draw_volcano(axes[0], fc, corrected, 'up',
                         f'{disp_label} — {cond_name} vs Naïve\nIncreased')
            draw_volcano(axes[1], fc, corrected, 'down',
                         f'{disp_label} — {cond_name} vs Naïve\nDecreased')

            fig.text(0.5, 0.01,
                     f'|FC| > {FC_THRESH}  |  Corrected p > {CORR_THRESH}'
                     f'  |  Axes: FC ±7, p 0–20{aw_note}',
                     ha='center', fontsize=10, color='#666666', style='italic')

            plt.tight_layout(pad=2.0, rect=[0, 0.04, 1, 1])
            pdf.savefig(fig, dpi=180, bbox_inches='tight')
            plt.close(fig)

            n_up   = int(((corrected > CORR_THRESH) & (fc >  FC_THRESH)).sum())
            n_down = int(((corrected > CORR_THRESH) & (fc < -FC_THRESH)).sum())
            print(f'{disp_label} — {cond_name}: ↑{n_up}  ↓{n_down}')

print(f'\nSaved: {OUT}')
