"""
Membrane fraction bar/dot plots — one panel per condition
==========================================================
Panel 1 (Intox): Slc8a2  — Naive vs Intoxication
Panel 2 (AW):    Sirt3   — Naive vs Acute Withdrawal
Panel 3 (PA):    Snca    — Naive vs Protracted Abstinence

Same purple pastel style as make_rep_dotplots.py.
AW-M-3 excluded from AW replicate columns.

Outputs:
  Membrane_Condition_BarPlots.pdf / .png
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'Membrane_Condition_BarPlots.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Membrane_Condition_BarPlots.png')

# ── Panel config ──────────────────────────────────────────────────────────────
# cond_prefix: column infix used to find replicate columns (e.g. 'I' → Memb_I-)
# corr_col:    corrected p-value column name in the sheet
PANELS = [
    dict(gene='Cat',    cond='Intox', cond_prefix='I',
         corr_col='Corrected',
         xlabel='Intox',
         title='Cat / Catalase\n(Membrane — Intox)'),
    dict(gene='Agap3',  cond='Intox', cond_prefix='I',
         corr_col='Corrected',
         xlabel='Intox',
         title='Agap3\n(Membrane — Intox)'),
    dict(gene='Sirt3',  cond='AW',    cond_prefix='AW',
         corr_col='Corrected.1',
         xlabel='Acute\nWithdrawal',
         title='Sirt3\n(Membrane — AW)'),
    dict(gene='Dmgdh',  cond='PA',    cond_prefix='PA',
         corr_col='Corrected.2',
         xlabel='Protracted\nAbstinence',
         title='Dmgdh\n(Membrane — PA)'),
]

SHEET  = 'Membrane'
PREFIX = 'Memb'

# ── Colors (same as rep dotplots) ─────────────────────────────────────────────
C_NAIVE = '#C9BBE8'
C_COND  = '#7B5EA7'
C_EDGE  = '#3D2B6B'

JITTER = 0.09
DOT_S  = 35
BAR_W  = 0.42

def asterisks(corr):
    if pd.isna(corr):      return 'ns'
    if corr >= 10:         return '****'
    if corr >= 7:          return '***'
    if corr >= 5:          return '**'
    if corr >= 3.3:        return '*'
    return 'ns'

# ── Load sheet once ───────────────────────────────────────────────────────────
np.random.seed(42)
df = pd.read_excel(DATA_FILE, sheet_name=SHEET)
df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

panels = []
for cfg in PANELS:
    row = df[df['Gene symbol'].str.lower() == cfg['gene'].lower()].iloc[0]

    n_cols   = [c for c in df.columns if f'{PREFIX}_N-' in str(c)]
    cond_cols = [c for c in df.columns
                 if f'{PREFIX}_{cfg["cond_prefix"]}-' in str(c)
                 and 'M-3' not in str(c)]   # exclude AW-M-3 for all (no-op for Intox/PA)

    naive = pd.to_numeric(row[n_cols],   errors='coerce').dropna().values
    cond  = pd.to_numeric(row[cond_cols], errors='coerce').dropna().values
    corr  = pd.to_numeric(row.get(cfg['corr_col'], np.nan), errors='coerce')

    panels.append(dict(
        title=cfg['title'],
        xlabel=cfg['xlabel'],
        naive=naive,
        cond=cond,
        corr=corr,
        ast=asterisks(corr),
    ))
    print(f"{cfg['gene']} ({cfg['cond']}): n_naive={len(naive)}, n_cond={len(cond)}, corr={corr:.2f}")

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
fig.patch.set_facecolor('white')

for ax, p in zip(axes, panels):
    naive, cond = p['naive'], p['cond']

    n_mean, n_sem = naive.mean(), stats.sem(naive)
    c_mean, c_sem = cond.mean(),  stats.sem(cond)

    # Filled bars
    ax.bar(0, n_mean, width=BAR_W, color=C_NAIVE,
           edgecolor=C_EDGE, linewidth=1.0, zorder=2)
    ax.bar(1, c_mean, width=BAR_W, color=C_COND,
           edgecolor=C_EDGE, linewidth=1.0, zorder=2)

    # SEM error bars
    ax.errorbar(0, n_mean, yerr=n_sem, fmt='none',
                color=C_EDGE, capsize=5, capthick=1.4, elinewidth=1.4, zorder=3)
    ax.errorbar(1, c_mean, yerr=c_sem, fmt='none',
                color=C_EDGE, capsize=5, capthick=1.4, elinewidth=1.4, zorder=3)

    # Individual dots
    jn = np.random.uniform(-JITTER, JITTER, len(naive))
    jc = np.random.uniform(-JITTER, JITTER, len(cond))
    ax.scatter(jn,     naive, s=DOT_S, color='white', edgecolors=C_EDGE,
               linewidths=0.8, zorder=4)
    ax.scatter(1 + jc, cond,  s=DOT_S, color='white', edgecolors=C_EDGE,
               linewidths=0.8, zorder=4)

    # Significance bracket
    all_vals = np.concatenate([naive, cond])
    y_range  = all_vals.max() - all_vals.min()
    bar_top  = max(0, n_mean, c_mean, all_vals.max())
    y_line   = bar_top + y_range * 0.12
    y_tick   = y_range * 0.04
    y_ast    = y_line  + y_range * 0.04

    ax.plot([0, 0, 1, 1],
            [y_line - y_tick, y_line, y_line, y_line - y_tick],
            color='#333333', lw=1.2, zorder=5)

    ast_str = p['ast']
    fs = 14 if ast_str != 'ns' else 11
    ax.text(0.5, y_ast, ast_str,
            ha='center', va='bottom', fontsize=fs, fontfamily='Arial',
            fontweight='bold' if ast_str != 'ns' else 'normal',
            color='#111111', zorder=5)

    # Axes styling
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Naïve', p['xlabel']], fontsize=13, fontfamily='Arial')
    if ax is axes[0]:
        ax.set_ylabel('Log2 Abundance', fontsize=15, fontfamily='Arial', labelpad=6)
    else:
        ax.set_ylabel('')
    ax.set_title(p['title'], fontsize=13, fontweight='bold',
                 fontfamily='Arial', pad=8, fontstyle='italic')
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.tick_params(axis='y', labelsize=13)
    ax.tick_params(axis='x', length=0, pad=6)

    cur_top = ax.get_ylim()[1]
    ax.set_ylim(top=max(cur_top, y_ast + y_range * 0.15))

fig.text(0.5, 0.01,
         '* corr. p > 3.3  |  ** > 5  |  *** > 7  |  **** > 10  '
         '(−log2 BH-corrected)  |  Bars = mean ± SEM  |  AW-M-3 excluded',
         ha='center', va='bottom', fontsize=7.5,
         color='#888888', style='italic', fontfamily='Arial')

plt.tight_layout(rect=[0, 0.04, 1, 1])

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
