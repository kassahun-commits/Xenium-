"""
Cytosol fraction box-and-whisker plots — 3 panels, all 4 conditions
====================================================================
Panel 1: Adh5  — Naïve / Intox / AW / PA  (Intox representative)
Panel 2: Apip  — Naïve / Intox / AW / PA  (AW representative)
Panel 3: Txn2  — Naïve / Intox / AW / PA  (PA representative)

Asterisks shown only for conditions with corrected p > 3.3 vs Naïve.
AW-M-3 excluded from AW replicates.

Outputs:
  Cytosol_BoxPlots.pdf / .png
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
OUT_PDF    = os.path.join(SCRIPT_DIR, 'Cytosol_BoxPlots.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Cytosol_BoxPlots.png')

SHEET  = 'Cytosol'
PREFIX = 'Cyto'

PANELS = [
    dict(gene='Adh5', title='Adh5\n(Cytosol)'),
    dict(gene='Apip', title='Apip\n(Cytosol)'),
    dict(gene='Txn2', title='Txn2\n(Cytosol)'),
]

CONDITIONS = [
    dict(label='Naïve', col_prefix='N',  corr_col=None,          color='#C9BBE8', exclude=None),
    dict(label='Intox', col_prefix='I',  corr_col='Corrected',   color='#A08CC0', exclude=None),
    dict(label='AW',    col_prefix='AW', corr_col='Corrected.1', color='#7B5EA7', exclude='M-3'),
    dict(label='PA',    col_prefix='PA', corr_col='Corrected.2', color='#5B3E87', exclude=None),
]

C_EDGE = '#3D2B6B'
JITTER = 0.09
DOT_S  = 35

def asterisks(corr):
    if pd.isna(corr):   return 'ns'
    if corr >= 10:      return '****'
    if corr >= 7:       return '***'
    if corr >= 5:       return '**'
    if corr >= 3.3:     return '*'
    return 'ns'

np.random.seed(42)
df = pd.read_excel(DATA_FILE, sheet_name=SHEET)
df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
fig.patch.set_facecolor('white')

for ax, cfg in zip(axes, PANELS):
    row = df[df['Gene symbol'].str.lower() == cfg['gene'].lower()].iloc[0]

    group_vals = []
    for cond in CONDITIONS:
        cols = [c for c in df.columns
                if f'{PREFIX}_{cond["col_prefix"]}-' in str(c)
                and (cond['exclude'] is None or cond['exclude'] not in str(c))]
        vals = pd.to_numeric(row[cols], errors='coerce').dropna().values
        group_vals.append(vals)
        print(f"{cfg['gene']} {cond['label']}: n={len(vals)}")

    positions = list(range(4))

    bp = ax.boxplot(group_vals, positions=positions, patch_artist=True,
                    widths=0.45, showfliers=False,
                    medianprops=dict(color=C_EDGE, linewidth=2),
                    whiskerprops=dict(color=C_EDGE, linewidth=1.2),
                    capprops=dict(color=C_EDGE, linewidth=1.2),
                    boxprops=dict(linewidth=1.2),
                    zorder=2)

    for patch, cond in zip(bp['boxes'], CONDITIONS):
        patch.set_facecolor(cond['color'])
        patch.set_edgecolor(C_EDGE)

    for i, (vals, cond) in enumerate(zip(group_vals, CONDITIONS)):
        jitter = np.random.uniform(-JITTER, JITTER, len(vals))
        ax.scatter(i + jitter, vals, s=DOT_S, color='white',
                   edgecolors=C_EDGE, linewidths=0.8, zorder=4)

    all_vals = np.concatenate(group_vals)
    y_range  = all_vals.max() - all_vals.min()
    y_max_per_group = [v.max() if len(v) > 0 else all_vals.max() for v in group_vals]

    for i, cond in enumerate(CONDITIONS[1:], 1):
        corr = pd.to_numeric(row.get(cond['corr_col'], np.nan), errors='coerce')
        ast  = asterisks(corr)
        if ast != 'ns':
            y_top = y_max_per_group[i] + y_range * 0.07
            ax.text(i, y_top, ast, ha='center', va='bottom',
                    fontsize=13, fontfamily='Arial', fontweight='bold',
                    color='#111111', zorder=5)

    ax.set_xlim(-0.6, 3.6)
    ax.set_xticks(positions)
    ax.set_xticklabels([c['label'] for c in CONDITIONS], fontsize=11, fontfamily='Arial')
    ax.set_ylabel('Log2 Abundance' if ax is axes[0] else '',
                  fontsize=13, fontfamily='Arial', labelpad=6)
    ax.set_title(cfg['title'], fontsize=12, fontweight='bold',
                 fontfamily='Arial', pad=8, fontstyle='italic')
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.tick_params(axis='y', labelsize=11)
    ax.tick_params(axis='x', length=0, pad=6)
    ax.set_ylim(top=max(ax.get_ylim()[1], max(y_max_per_group) + y_range * 0.28))

fig.text(0.5, 0.01,
         '* corr. p > 3.3  |  ** > 5  |  *** > 7  |  **** > 10  '
         '(−log2 BH-corrected)  |  Box = median ± IQR  |  AW-M-3 excluded',
         ha='center', va='bottom', fontsize=7.5,
         color='#888888', style='italic', fontfamily='Arial')

plt.tight_layout(rect=[0, 0.04, 1, 1])

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
