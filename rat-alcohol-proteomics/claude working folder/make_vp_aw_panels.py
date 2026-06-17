import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'VP_AcuteWithdrawal_UpDown.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_UP    = '#E8305A'   # bold rose/magenta
C_DOWN  = '#2B7FD4'   # bold blue
C_OTHER = '#C8C8C8'   # greyed-out significant (wrong direction)
C_NS    = '#EBEBEB'   # non-significant

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_N-',  'Memb_AW-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-',  'Cyto_AW-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_N-', 'Chrom_AW-', 'keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_N-',   'Nuc_AW-',   'keep_review'),
]

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        if len(n) >= 2 and len(c) >= 2:
            _, p = stats.ttest_ind(n, c, equal_var=False)
            p_vals.append(p)
        else:
            p_vals.append(np.nan)
    p_vals    = pd.Series(p_vals, index=df.index)
    valid     = p_vals.notna()
    ranks     = p_vals[valid].rank(ascending=False)
    corrected = pd.Series(np.nan, index=df.index)
    corrected[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corrected

def make_panel(ax, fc, corrected, highlight, title, c_highlight, c_other):
    """highlight: 'up' or 'down'"""
    mask = fc.notna() & corrected.notna()
    f  = fc[mask]
    cr = corrected[mask]

    up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
    down = (cr > CORR_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down

    if highlight == 'up':
        hi  = up
        lo  = down
    else:
        hi  = down
        lo  = up

    # Plot order: NS first, then wrong-direction sig, then highlighted on top
    ax.scatter(f[ns],  cr[ns],  s=18,  color=C_NS,    alpha=0.5,  linewidths=0, rasterized=True)
    ax.scatter(f[lo],  cr[lo],  s=30,  color=C_OTHER, alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(f[hi],  cr[hi],  s=60,  color=c_highlight, alpha=1.0, linewidths=0, rasterized=True)

    # Threshold lines
    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlabel('Fold Change', fontsize=20, labelpad=10)
    ax.set_ylabel('Corrected p-value\n(−log₂ adjusted p)', fontsize=20, labelpad=10)
    ax.set_title(title, fontsize=22, fontweight='bold', pad=12)
    ax.tick_params(labelsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    n_hi  = int(hi.sum())
    n_tot = int(mask.sum())
    direction = '↑' if highlight == 'up' else '↓'
    ax.text(0.98, 0.98, f'{direction}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=18, color='#333333', fontweight='bold')

    dir_label = 'Increased' if highlight == 'up' else 'Decreased'
    ax.legend(handles=[
        mpatches.Patch(color=c_highlight, label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_OTHER,     label='Other sig.'),
        mpatches.Patch(color=C_NS,        label='NS'),
    ], fontsize=14, loc='upper left', framealpha=0.85,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.9)

with PdfPages(OUT) as pdf:
    for label, excel_name, naive_pat, aw_pat, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
        if fmode == 'keep_review':
            df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        naive_cols = [c for c in df.columns if naive_pat in str(c)]
        aw_cols    = [c for c in df.columns if aw_pat    in str(c)]
        fc, corrected = calc_stats(df, naive_cols, aw_cols)

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        make_panel(axes[0], fc, corrected, 'up',
                   f'{label}\nAcute Withdrawal vs Naïve — Increased',
                   C_UP, C_DOWN)

        make_panel(axes[1], fc, corrected, 'down',
                   f'{label}\nAcute Withdrawal vs Naïve — Decreased',
                   C_DOWN, C_UP)

        fig.text(0.5, 0.01,
                 f'Thresholds: |Fold Change| > {FC_THRESH}  |  Corrected p-value > {CORR_THRESH}',
                 ha='center', fontsize=13, color='#555555', style='italic')

        plt.tight_layout(pad=2.0, rect=[0, 0.04, 1, 1])
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        n_up   = int(((corrected > CORR_THRESH) & (fc >  FC_THRESH)).sum())
        n_down = int(((corrected > CORR_THRESH) & (fc < -FC_THRESH)).sum())
        print(f'{label}: ↑{n_up}  ↓{n_down}')

print(f'Saved: {OUT}')
