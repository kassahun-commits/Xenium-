"""
4-panel dot/bar plots: Naive vs Intox rep values
Genes: Slc8a1 (Membrane), Scn1b (Membrane), Rack1 (Chromatin), Adh5 (Cytosol)
Y-axis: log2 abundance (rep values as stored in sheet)
Significance asterisks from pre-computed corrected p-value (same threshold as project)
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
OUT_PDF    = os.path.join(SCRIPT_DIR, 'CandidateGenes_RepBarPlots_4panel.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'CandidateGenes_RepBarPlots_4panel.png')

# ── Gene config — 4 panel version (Membrane, Cytosol, SN, Chromatin) ──────────
GENES = [
    dict(gene='Slc8a1',  sheet='Membrane',        prefix='Memb',  filt=False,
         title='Slc8a1 / NCX1\n(Membrane)'),
    dict(gene='Adh5',    sheet='Cytosol',         prefix='Cyto',  filt=False,
         title='Adh5\n(Cytosol)'),
    dict(gene='Ndufa13', sheet='Soluble nuclear', prefix='Nuc',   filt=False,
         title='Ndufa13 / GRIM-19\n(Soluble Nuclear)'),
    dict(gene='Rack1',   sheet='Chromatin',       prefix='Chrom', filt=True,
         title='Rack1\n(Chromatin)'),
]

# Colors — pastel purple palette
C_NAIVE = '#C9BBE8'   # light lavender
C_INTOX = '#7B5EA7'   # medium purple
C_EDGE  = '#3D2B6B'   # dark purple edge

JITTER  = 0.09
DOT_S   = 35          # smaller dots
BAR_W   = 0.42
CAP_W   = 0.15

def asterisks(corr):
    """Convert corrected p (−log2 scale) to asterisk string."""
    if corr >= 10: return '****'
    if corr >= 7:  return '***'
    if corr >= 5:  return '**'
    if corr >= 3.3: return '*'
    return 'ns'

# ── Load data ──────────────────────────────────────────────────────────────────
np.random.seed(42)
sheet_cache = {}

def get_sheet(sheet, filt):
    key = (sheet, filt)
    if key not in sheet_cache:
        df = pd.read_excel(DATA_FILE, sheet_name=sheet)
        df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
        if filt and 'Filter' in df.columns:
            df = df[df['Filter'].isin(['Keep', 'Review'])]
        sheet_cache[key] = df
    return sheet_cache[key]

panels = []
for cfg in GENES:
    df  = get_sheet(cfg['sheet'], cfg['filt'])
    row = df[df['Gene symbol'].str.lower() == cfg['gene'].lower()].iloc[0]

    n_cols = [c for c in df.columns if f"{cfg['prefix']}_N-" in str(c)]
    i_cols = [c for c in df.columns if f"{cfg['prefix']}_I-" in str(c)]

    naive = pd.to_numeric(row[n_cols], errors='coerce').dropna().values
    intox = pd.to_numeric(row[i_cols], errors='coerce').dropna().values
    corr  = pd.to_numeric(row['Corrected'], errors='coerce')

    panels.append(dict(
        title=cfg['title'],
        naive=naive, intox=intox,
        corr=corr,
        ast=asterisks(corr),
    ))

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
fig.patch.set_facecolor('white')

for ax, p in zip(axes, panels):
    naive, intox = p['naive'], p['intox']

    # Means and SEM
    n_mean, n_sem = naive.mean(), stats.sem(naive)
    i_mean, i_sem = intox.mean(), stats.sem(intox)

    # ── Filled bars (0 → mean) ──
    ax.bar(0, n_mean, width=BAR_W, color=C_NAIVE,
           edgecolor=C_EDGE, linewidth=1.0, zorder=2)
    ax.bar(1, i_mean, width=BAR_W, color=C_INTOX,
           edgecolor=C_EDGE, linewidth=1.0, zorder=2)

    # ── SEM error bars ──
    ax.errorbar(0, n_mean, yerr=n_sem, fmt='none',
                color=C_EDGE, capsize=5, capthick=1.4,
                elinewidth=1.4, zorder=3)
    ax.errorbar(1, i_mean, yerr=i_sem, fmt='none',
                color=C_EDGE, capsize=5, capthick=1.4,
                elinewidth=1.4, zorder=3)

    # ── Individual dots on top ──
    jn = np.random.uniform(-JITTER, JITTER, len(naive))
    ji = np.random.uniform(-JITTER, JITTER, len(intox))
    ax.scatter(jn,     naive, s=DOT_S, color='white', edgecolors=C_EDGE,
               linewidths=0.8, zorder=4)
    ax.scatter(1 + ji, intox, s=DOT_S, color='white', edgecolors=C_EDGE,
               linewidths=0.8, zorder=4)

    # ── Significance bracket ──
    all_vals = np.concatenate([naive, intox])
    y_range  = all_vals.max() - all_vals.min()
    # Bracket must clear the top of both bars (bars go from 0, so top = max(0, means, data))
    bar_top  = max(0, n_mean, i_mean, all_vals.max())
    y_line   = bar_top + y_range * 0.12
    y_tick   = y_range * 0.04
    y_ast    = y_line + y_range * 0.04

    ax.plot([0, 0, 1, 1],
            [y_line - y_tick, y_line, y_line, y_line - y_tick],
            color='#333333', lw=1.2, zorder=5)

    ast_str = p['ast']
    fs = 14 if ast_str != 'ns' else 11
    ax.text(0.5, y_ast, ast_str,
            ha='center', va='bottom',
            fontsize=fs, fontfamily='Arial',
            fontweight='bold' if ast_str != 'ns' else 'normal',
            color='#111111', zorder=5)

    # ── Axes styling ──
    ax.set_xlim(-0.6, 1.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Naïve', 'Intox'], fontsize=15, fontfamily='Arial')
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

    # Extend y-axis to fit bracket + asterisks
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
