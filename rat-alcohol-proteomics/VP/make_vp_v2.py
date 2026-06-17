import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE = 'EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Volcano_Plots_v2.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# Pastel colour scheme
C_UP   = '#F4A0A0'   # soft rose
C_DOWN = '#A0BCF4'   # soft periwinkle
C_NS   = '#CCCCCC'   # light grey

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
    p_vals   = pd.Series(p_vals, index=df.index)
    valid    = p_vals.notna()
    ranks    = p_vals[valid].rank(ascending=False)
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

    ax.scatter(cr[ns],   f[ns],   s=18,  color=C_NS,   alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(cr[down], f[down], s=40,  color=C_DOWN, alpha=0.85, linewidths=0, rasterized=True)
    ax.scatter(cr[up],   f[up],   s=40,  color=C_UP,   alpha=0.85, linewidths=0, rasterized=True)

    ax.axvline(CORR_THRESH, color='#666666', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.axhline( FC_THRESH,  color='#666666', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.axhline(-FC_THRESH,  color='#666666', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.6, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.6, alpha=0.4)

    ax.set_xlabel('Corrected p-value (−log₂ adjusted p)', fontsize=22, labelpad=10)
    ax.set_ylabel('Fold Change', fontsize=22, labelpad=10)
    ax.set_title(title, fontsize=22, fontweight='bold', pad=12)
    ax.tick_params(labelsize=17)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)

    n_up   = int(up.sum())
    n_down = int(down.sum())
    n_tot  = int(mask.sum())
    ax.text(0.98, 0.98, f'↑{n_up}   ↓{n_down}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top', fontsize=17, color='#444444')
    ax.legend(handles=[
        mpatches.Patch(color=C_UP,   label=f'Up  (n={n_up})'),
        mpatches.Patch(color=C_DOWN, label=f'Down  (n={n_down})'),
        mpatches.Patch(color=C_NS,   label='NS'),
    ], fontsize=15, loc='lower right', framealpha=0.8, edgecolor='#CCCCCC',
       frameon=True, borderpad=0.8)

CONDITIONS = [
    ('Intoxication vs Naïve',          'I'),
    ('Acute Withdrawal vs Naïve',      'AW'),
    ('Protracted Abstinence vs Naïve', 'PA'),
]

SHEETS = [
    ('Chromatin',       'Chrom_N-', 'Chrom_{cond}-'),
    ('Soluble nuclear', 'Nuc_N-',   'Nuc_{cond}-'),
]

# Build job list: (title, sheet_name, filter_mode, naive_pat, cond_pat)
jobs = []
for sheet_name, naive_pat, cond_pat_tmpl in SHEETS:
    for filter_label, filter_mode in [('Keep only', 'keep'), ('Keep + Review', 'keep_review')]:
        for cond_label, cond_code in CONDITIONS:
            cond_pat = cond_pat_tmpl.replace('{cond}', cond_code)
            title = f'{sheet_name}\n{cond_label}\n({filter_label})'
            jobs.append((title, sheet_name, filter_mode, naive_pat, cond_pat))

# Cache sheet data
sheet_cache = {}
for _, sheet_name, _, _, _ in jobs:
    if sheet_name not in sheet_cache:
        sheet_cache[sheet_name] = pd.read_excel(FILE, sheet_name=sheet_name, header=0)

with PdfPages(OUT) as pdf:
    for page_start in range(0, len(jobs), 3):
        page_jobs = jobs[page_start:page_start + 3]
        n_plots = len(page_jobs)
        fig, axes = plt.subplots(1, n_plots, figsize=(10 * n_plots, 9))
        if n_plots == 1:
            axes = [axes]

        for ax, (title, sheet_name, filter_mode, naive_pat, cond_pat) in zip(axes, page_jobs):
            df = sheet_cache[sheet_name].copy()
            if filter_mode == 'keep':
                df = df[df['Filter'] == 'Keep'].reset_index(drop=True)
            elif filter_mode == 'keep_review':
                df = df[df['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)

            naive_cols = [c for c in df.columns if naive_pat in str(c)]
            cond_cols  = [c for c in df.columns if cond_pat  in str(c)]
            fc, corrected = calc_stats(df, naive_cols, cond_cols)
            make_plot(ax, fc, corrected, title)

        fig.text(0.5, 0.01,
                 f'Thresholds: Corrected p-value > {CORR_THRESH} (vertical dashed)  |  '
                 f'Fold Change > ±{FC_THRESH} (horizontal dashed)',
                 ha='center', fontsize=16, color='#555555', style='italic')
        plt.tight_layout(pad=2.5, rect=[0, 0.04, 1, 1])
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        page_num = page_start // 3 + 1
        print(f'Page {page_num}: {[j[0].replace(chr(10), " | ") for j in page_jobs]}')

print(f'\nSaved: {OUT}')
print(f'Total plots: {len(jobs)}')
