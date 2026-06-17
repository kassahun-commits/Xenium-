import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Volcano_Plots_v4.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_UP   = '#E8735A'
C_DOWN = '#6B9E78'
C_NS   = '#C8C8C8'

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

def make_plot(ax, fc, corrected, title):
    mask = fc.notna() & corrected.notna()
    f  = fc[mask]
    cr = corrected[mask]

    up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
    down = (cr > CORR_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down

    # X = Fold Change, Y = Corrected p-value
    ax.scatter(f[ns],   cr[ns],   s=30,  color=C_NS,   alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(f[down], cr[down], s=60,  color=C_DOWN, alpha=0.85, linewidths=0, rasterized=True)
    ax.scatter(f[up],   cr[up],   s=60,  color=C_UP,   alpha=0.85, linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#666666', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#666666', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#666666', linestyle='--', linewidth=2.0, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=1.0, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=1.0, alpha=0.4)

    ax.set_xlabel('Fold Change', fontsize=56, labelpad=16)
    ax.set_ylabel('Corrected p-value\n(−log₂ adjusted p)', fontsize=56, labelpad=16)
    ax.set_title(title, fontsize=58, fontweight='bold', pad=20)
    ax.tick_params(labelsize=48)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2.5)
    ax.spines['bottom'].set_linewidth(2.5)

    n_up   = int(up.sum())
    n_down = int(down.sum())
    n_tot  = int(mask.sum())
    ax.text(0.98, 0.98, f'n={n_tot}',
            transform=ax.transAxes, ha='right', va='top', fontsize=52, color='#444444',
            fontweight='bold')
    ax.legend(handles=[
        mpatches.Patch(color=C_UP,   label=f'Up  (n={n_up})'),
        mpatches.Patch(color=C_DOWN, label=f'Down  (n={n_down})'),
        mpatches.Patch(color=C_NS,   label='NS'),
    ], fontsize=44, loc='upper left', framealpha=0.8, edgecolor='#CCCCCC', frameon=True, borderpad=0.8)

CONDITIONS = [
    ('Intoxication vs Naïve',          'I'),
    ('Acute Withdrawal vs Naïve',      'AW'),
    ('Protracted Abstinence vs Naïve', 'PA'),
]

# sheet_name, naive_pat, cond_pat_tmpl, filter_mode ('keep_review' or 'all')
# Order: Membrane, Cytosol, Chromatin, Soluble nuclear
SHEETS = [
    ('Membrane',        'Memb_N-',  'Memb_{cond}-',  'all'),
    ('Cytosol',         'Cyto_N-',  'Cyto_{cond}-',  'all'),
    ('Chromatin',       'Chrom_N-', 'Chrom_{cond}-', 'keep_review'),
    ('Soluble nuclear', 'Nuc_N-',   'Nuc_{cond}-',   'keep_review'),
]

sheet_cache = {s[0]: pd.read_excel(FILE, sheet_name=s[0], header=0) for s in SHEETS}

# ── PDF: 4 compartments × 3 conditions = 12 plots, 4 pages ──
jobs = []
for sheet_name, naive_pat, cond_pat_tmpl, fmode in SHEETS:
    for cond_label, cond_code in CONDITIONS:
        cond_pat = cond_pat_tmpl.replace('{cond}', cond_code)
        title = f'{sheet_name}\n{cond_label}'
        jobs.append((title, sheet_name, fmode, naive_pat, cond_pat))

for out_file in [OUT, 'Volcano_Plots_Alt_Colors.pdf']:
    with PdfPages(out_file) as pdf:
        for page_start in range(0, len(jobs), 3):
            page_jobs = jobs[page_start:page_start + 3]
            fig, axes = plt.subplots(1, 3, figsize=(54, 20))
            for ax, (title, sheet_name, fmode, naive_pat, cond_pat) in zip(axes, page_jobs):
                df = sheet_cache[sheet_name].copy()
                if fmode == 'keep_review':
                    df = df[df['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
                naive_cols = [c for c in df.columns if naive_pat in str(c)]
                cond_cols  = [c for c in df.columns if cond_pat  in str(c)]
                fc, corrected = calc_stats(df, naive_cols, cond_cols)
                make_plot(ax, fc, corrected, title)
            fig.text(0.5, 0.01,
                     f'Thresholds: |Fold Change| > {FC_THRESH} (vertical dashed)  |  '
                     f'Corrected p-value > {CORR_THRESH} (horizontal dashed)',
                     ha='center', fontsize=40, color='#555555', style='italic')
            plt.tight_layout(pad=2.5, rect=[0, 0.04, 1, 1])
            pdf.savefig(fig, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f'{out_file} — page {page_start//3+1} done')
    print(f'Saved: {out_file}')
