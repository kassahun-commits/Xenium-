import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE = 'EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Volcano_Plots_FC_vs_Score.pdf'

SCORE_THRESH = 1.65   # x-axis threshold
FC_THRESH    = 0.5    # y-axis threshold (fold change)

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        if len(n) >= 2 and len(c) >= 2:
            _, p = stats.ttest_ind(n, c, equal_var=False)
            p_vals.append(p)
        else:
            p_vals.append(np.nan)
    p_vals = pd.Series(p_vals, index=df.index)
    valid = p_vals.notna()
    ranks_desc = p_vals[valid].rank(ascending=False)
    corrected = pd.Series(np.nan, index=df.index)
    corrected[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks_desc))
    score = fc * corrected
    return fc, score

def make_plot(ax, fc, score, title):
    mask = fc.notna() & score.notna()
    f   = fc[mask]
    sc  = score[mask]

    up   = (sc >  SCORE_THRESH) & (f >  FC_THRESH)
    down = (sc < -SCORE_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down

    ax.scatter(sc[ns],   f[ns],   s=6,  color='#AAAAAA', alpha=0.4, linewidths=0, rasterized=True)
    ax.scatter(sc[down], f[down], s=12, color='#2C7BB6', alpha=0.85, linewidths=0, rasterized=True)
    ax.scatter(sc[up],   f[up],   s=12, color='#D7191C', alpha=0.85, linewidths=0, rasterized=True)

    ax.axvline( SCORE_THRESH, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axvline(-SCORE_THRESH, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline( FC_THRESH,    color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(-FC_THRESH,    color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(0,             color='black', linestyle='-',  linewidth=0.4, alpha=0.3)
    ax.axvline(0,             color='black', linestyle='-',  linewidth=0.4, alpha=0.3)

    ax.set_xlabel('Score (Fold Change × Corrected)', fontsize=9)
    ax.set_ylabel('Fold Change', fontsize=9)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=6)
    ax.tick_params(labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_up = up.sum(); n_down = down.sum()
    ax.text(0.98, 0.98, f'↑{n_up}  ↓{n_down}  n={mask.sum()}',
            transform=ax.transAxes, ha='right', va='top', fontsize=7, color='#333333')
    ax.legend(handles=[
        mpatches.Patch(color='#D7191C', label=f'Up ({n_up})'),
        mpatches.Patch(color='#2C7BB6', label=f'Down ({n_down})'),
        mpatches.Patch(color='#AAAAAA', label='NS'),
    ], fontsize=6, loc='lower right', framealpha=0.7, edgecolor='none')

CONDITIONS = [
    ('Intoxication vs Naïve',          '_I-'),
    ('Acute Withdrawal vs Naïve',      '_AW-'),
    ('Protracted Abstinence vs Naïve', '_PA-'),
]

df_all = pd.read_excel(FILE, sheet_name='Chromatin', header=0)
naive_cols = [c for c in df_all.columns if 'Chrom_N-' in str(c)]

jobs = []
for filter_label, filter_mode in [('Keep only', 'keep'), ('Keep + Review', 'keep_review')]:
    for cond_label, cond_pat in CONDITIONS:
        cond_cols = [c for c in df_all.columns if cond_pat in str(c) and 'Chrom' in str(c)]
        jobs.append((f'Chromatin — {cond_label}\n({filter_label})', filter_mode, naive_cols, cond_cols))

with PdfPages(OUT) as pdf:
    for page_start in range(0, len(jobs), 3):
        page_jobs = jobs[page_start:page_start+3]
        fig, axes = plt.subplots(1, 3, figsize=(16.5, 5))
        for ax, (title, filter_mode, naive_cols, cond_cols) in zip(axes, page_jobs):
            df = df_all.copy()
            if filter_mode == 'keep':
                df = df[df['Filter'] == 'Keep'].reset_index(drop=True)
            elif filter_mode == 'keep_review':
                df = df[df['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
            fc, score = calc_stats(df, naive_cols, cond_cols)
            make_plot(ax, fc, score, title)

        # Add threshold note at bottom of page
        fig.text(0.5, 0.01,
                 f'Thresholds: Score > ±{SCORE_THRESH} (vertical dashed) | Fold Change > ±{FC_THRESH} (horizontal dashed) — Score = Fold Change × Corrected',
                 ha='center', fontsize=7, color='#555555', style='italic')
        plt.tight_layout(pad=2.0, rect=[0, 0.03, 1, 1])
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'Page {page_start//3 + 1} done')

print(f'\nSaved: {OUT}')
