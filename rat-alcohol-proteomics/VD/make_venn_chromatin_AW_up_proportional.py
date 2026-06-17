"""
Proportional Venn diagram — Chromatin AW Upregulated proteins
=============================================================
AW Increased vs Intox Increased vs PA Increased (Chromatin fraction).
Circle areas proportional to actual set sizes.
AW-M-3 excluded. Chromatin: Keep+Review only.

Output: Venn_Chromatin_AW_Up_Proportional.pdf / .png
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE    = os.path.join(SCRIPT_DIR,
          '../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_PDF = os.path.join(SCRIPT_DIR, 'Venn_Chromatin_AW_Up_Proportional.pdf')
OUT_PNG = os.path.join(SCRIPT_DIR, 'Venn_Chromatin_AW_Up_Proportional.png')

CORR_THRESH = 3.3
FC_THRESH   = 0.5

N_SUFFIXES  = ['N-F-1',  'N-F-2',  'N-F-3',  'N-M-1',  'N-M-2']
I_SUFFIXES  = ['I-F-1',  'I-F-2',  'I-F-3',  'I-M-1',  'I-M-2']
AW_SUFFIXES = ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2']   # AW-M-3 excluded
PA_SUFFIXES = ['PA-F-1', 'PA-F-2', 'PA-F-3', 'PA-M-1', 'PA-M-2']

PREFIX = 'Chrom_'

# ── Load Chromatin ─────────────────────────────────────────────────────────────
raw = pd.read_excel(FILE, sheet_name='Chromatin')
df  = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

n_cols  = [c for c in df.columns if any(PREFIX+s in str(c) for s in N_SUFFIXES)]
i_cols  = [c for c in df.columns if any(PREFIX+s in str(c) for s in I_SUFFIXES)]
aw_cols = [c for c in df.columns if any(PREFIX+s in str(c) for s in AW_SUFFIXES)]
pa_cols = [c for c in df.columns if any(PREFIX+s in str(c) for s in PA_SUFFIXES)]

# ── Significance ───────────────────────────────────────────────────────────────
def get_sig_genes(naive_cols, cond_cols, direction):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pvs   = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pvs.append(stats.ttest_ind(n, c, equal_var=False)[1]
                   if len(n) >= 2 and len(c) >= 2 else np.nan)
    pvs   = pd.Series(pvs, index=df.index)
    valid = pvs.notna()
    ranks = pvs[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pvs[valid] * valid.sum() / ranks))
    mask = (corr > CORR_THRESH) & (fc > FC_THRESH if direction == 'up' else fc < -FC_THRESH)
    return set(df.loc[mask, 'Gene symbol'].tolist())

AW    = get_sig_genes(n_cols, aw_cols, 'up')
INTOX = get_sig_genes(n_cols, i_cols,  'up')
PA    = get_sig_genes(n_cols, pa_cols, 'up')

# Actual subset sizes (for labels)
def actual_subsets(A, B, C):
    return (
        len(A - B - C),      # A only
        len(B - A - C),      # B only
        len((A & B) - C),    # A & B only
        len(C - A - B),      # C only
        len((A & C) - B),    # A & C only
        len((B & C) - A),    # B & C only
        len(A & B & C),      # all
    )

actual = actual_subsets(AW, INTOX, PA)
print(f'AW up:    n={len(AW)}')
print(f'Intox up: n={len(INTOX)}')
print(f'PA up:    n={len(PA)}')
print(f'Subsets (AW-only, Intox-only, AW∩Intox, PA-only, AW∩PA, Intox∩PA, all): {actual}')

# ── Colors ────────────────────────────────────────────────────────────────────
import matplotlib.colors as mc

C_AW    = '#5FAD8E'   # teal green
C_INTOX = '#D4845A'   # salmon orange
C_PA    = '#9B8EC4'   # lavender purple

def blend2(c1, c2):
    r1, r2 = np.array(mc.to_rgb(c1)), np.array(mc.to_rgb(c2))
    return mc.to_hex((r1 + r2) / 2)

def blend3(c1, c2, c3):
    r1, r2, r3 = np.array(mc.to_rgb(c1)), np.array(mc.to_rgb(c2)), np.array(mc.to_rgb(c3))
    return mc.to_hex((r1 + r2 + r3) / 3)

PATCH_COLORS = {
    '100': C_AW,
    '010': C_INTOX,
    '001': C_PA,
    '110': blend2(C_AW, C_INTOX),
    '101': blend2(C_AW, C_PA),
    '011': blend2(C_INTOX, C_PA),
    '111': blend3(C_AW, C_INTOX, C_PA),
}

# ── Layout: sqrt-scaled so circles form a balanced triangle ───────────────────
# Actual: (AW-only, Intox-only, AW∩Intox, PA-only, AW∩PA, Intox∩PA, all)
layout = tuple(max(3.0, v ** 0.5) for v in actual)

# ── Draw ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 9))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

v = venn3(subsets=layout,
          set_labels=(
              f'AW Increased\n(n={len(AW)})',
              f'Intox. Increased\n(n={len(INTOX)})',
              f'PA Increased\n(n={len(PA)})',
          ),
          ax=ax)

# Color each patch
for pid, color in PATCH_COLORS.items():
    patch = v.get_patch_by_id(pid)
    if patch:
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_edgecolor('none')

# Circle outlines
vc = venn3_circles(subsets=layout, ax=ax)
for circle in (vc.patches if hasattr(vc, 'patches') else vc):
    circle.set_edgecolor('#333333')
    circle.set_linewidth(2.5)
    circle.set_facecolor('none')

# Override labels with actual counts
for lbl, val in zip(v.subset_labels, actual):
    if lbl is not None:
        lbl.set_text(str(val))
        lbl.set_fontsize(21)
        lbl.set_fontfamily('Arial')
        lbl.set_fontweight('bold')
        lbl.set_color('#111111')

# Set labels
for lbl in v.set_labels:
    if lbl:
        lbl.set_fontsize(19)
        lbl.set_fontfamily('Arial')
        lbl.set_fontweight('bold')
        lbl.set_color('#111111')

ax.set_title('Chromatin — AW Upregulated Proteins\nOverlap Across Conditions',
             fontsize=18, fontweight='bold', fontfamily='Arial', pad=16)

fig.text(0.5, 0.02,
         f'|FC| > {FC_THRESH}, corrected p > {CORR_THRESH}  |  AW-M-3 excluded  |  Chromatin: Keep+Review only',
         ha='center', va='bottom', fontsize=9,
         color='#888888', style='italic', fontfamily='Arial')

plt.tight_layout(rect=[0, 0.04, 1, 1])

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
