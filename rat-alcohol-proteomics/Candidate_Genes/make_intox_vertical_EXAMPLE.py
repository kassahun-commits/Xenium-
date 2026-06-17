"""
EXAMPLE — single panel (Cat / Membrane) with enlarged fonts.
Shows what the vertical Intox figure will look like before committing.
Output: Intox_Vertical_EXAMPLE.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Intox_Vertical_EXAMPLE.png')

# ── Font sizes (triple the original ~11-13 pt) ────────────────────────────────
FS_TITLE   = 36    # panel title
FS_TICK    = 30    # x/y tick labels
FS_YLABEL  = 32    # y-axis label
FS_AST     = 38    # asterisk
FS_FOOT    = 20    # footnote

CONDITIONS = [
    dict(label='Naïve', col_prefix='N',  corr_col=None,          color='#C9BBE8', exclude=None),
    dict(label='Intox', col_prefix='I',  corr_col='Corrected',   color='#A08CC0', exclude=None),
    dict(label='AW',    col_prefix='AW', corr_col='Corrected.1', color='#7B5EA7', exclude='M-3'),
    dict(label='PA',    col_prefix='PA', corr_col='Corrected.2', color='#5B3E87', exclude=None),
]
C_EDGE = '#3D2B6B'
JITTER = 0.09
DOT_S  = 120    # bigger dots to match bigger figure

def asterisks(corr):
    if pd.isna(corr):  return 'ns'
    if corr >= 10:     return '****'
    if corr >= 7:      return '***'
    if corr >= 5:      return '**'
    if corr >= 3.3:    return '*'
    return 'ns'

np.random.seed(42)
df = pd.read_excel(DATA_FILE, sheet_name='Membrane')
df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
row = df[df['Gene symbol'].str.lower() == 'cat'].iloc[0]
PREFIX = 'Memb'

group_vals = []
for cond in CONDITIONS:
    cols = [c for c in df.columns
            if f'{PREFIX}_{cond["col_prefix"]}-' in str(c)
            and (cond['exclude'] is None or cond['exclude'] not in str(c))]
    vals = pd.to_numeric(row[cols], errors='coerce').dropna().values
    group_vals.append(vals)

fig, ax = plt.subplots(1, 1, figsize=(7, 6))
fig.patch.set_facecolor('white')

bp = ax.boxplot(group_vals, positions=list(range(4)), patch_artist=True,
                widths=0.45, showfliers=False,
                medianprops=dict(color=C_EDGE, linewidth=2.5),
                whiskerprops=dict(color=C_EDGE, linewidth=1.8),
                capprops=dict(color=C_EDGE, linewidth=1.8),
                boxprops=dict(linewidth=1.8), zorder=2)

for patch, cond in zip(bp['boxes'], CONDITIONS):
    patch.set_facecolor(cond['color'])
    patch.set_edgecolor(C_EDGE)

for i, (vals, cond) in enumerate(zip(group_vals, CONDITIONS)):
    jitter = np.random.uniform(-JITTER, JITTER, len(vals))
    ax.scatter(i + jitter, vals, s=DOT_S, color='white',
               edgecolors=C_EDGE, linewidths=1.2, zorder=4)

all_vals   = np.concatenate(group_vals)
y_range    = all_vals.max() - all_vals.min()
y_max_pg   = [v.max() if len(v) > 0 else all_vals.max() for v in group_vals]

for i, cond in enumerate(CONDITIONS[1:], 1):
    corr = pd.to_numeric(row.get(cond['corr_col'], np.nan), errors='coerce')
    ast  = asterisks(corr)
    if ast != 'ns':
        y_top = y_max_pg[i] + y_range * 0.07
        ax.text(i, y_top, ast, ha='center', va='bottom',
                fontsize=FS_AST, fontfamily='Arial', fontweight='bold',
                color='#111111', zorder=5)

ax.set_xlim(-0.6, 3.6)
ax.set_xticks(list(range(4)))
ax.set_xticklabels([c['label'] for c in CONDITIONS], fontsize=FS_TICK, fontfamily='Arial')
ax.set_ylabel('Log2 Abundance', fontsize=FS_YLABEL, fontfamily='Arial', labelpad=6)
ax.set_title('Cat / Catalase\n(Membrane)', fontsize=FS_TITLE, fontweight='bold',
             fontfamily='Arial', pad=10, fontstyle='italic')
ax.set_facecolor('white')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
ax.tick_params(axis='y', labelsize=FS_TICK)
ax.tick_params(axis='x', length=0, pad=8)
ax.set_ylim(top=max(ax.get_ylim()[1], max(y_max_pg) + y_range * 0.32))

plt.tight_layout()
fig.savefig(OUT_PNG, dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PNG}')
