"""
Venn diagrams — AW-M-3 excluded from Acute Withdrawal.
Produces:
  Venn_AcuteWithdrawal_noAWM3.pdf   — one page per compartment (4 pages)
  Venn_Summary_noAWM3.pdf           — all 8 Venns on one page (4 rows × 2 cols)

Filters:
  Chromatin       — Keep + Review only
  All others      — all rows (no filter)
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_INDIVIDUAL = 'Venn_AcuteWithdrawal_noAWM3.pdf'
OUT_SUMMARY    = 'Venn_Summary_noAWM3.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_UP   = '#F28FAD'
C_DOWN = '#7BAFD4'

# AW-M-3 explicitly excluded
AW_SUFFIXES = ['AW-F-1', 'AW-F-2', 'AW-M-1', 'AW-M-2']
I_SUFFIXES  = ['I-F-1',  'I-F-2',  'I-F-3',  'I-M-1',  'I-M-2']
PA_SUFFIXES = ['PA-F-1', 'PA-F-2', 'PA-F-3', 'PA-M-1', 'PA-M-2']
N_SUFFIXES  = ['N-F-1',  'N-F-2',  'N-F-3',  'N-M-1',  'N-M-2']

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'all'),
]

# ── Stats helpers ──────────────────────────────────────────────────────────────
def calc_stats(df, naive_cols, cond_cols):
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
    return fc, corr

def get_sig_genes(df, naive_cols, cond_cols, direction):
    fc, corr = calc_stats(df, naive_cols, cond_cols)
    mask = fc.notna() & corr.notna() & (corr > CORR_THRESH)
    sig  = mask & (fc > FC_THRESH if direction == 'up' else fc < -FC_THRESH)
    return set(df.loc[sig, 'Gene symbol'].astype(str).tolist())

def subset_sizes(A, B, C):
    return (
        len(A - B - C), len(B - A - C), len((A & B) - C),
        len(C - A - B), len((A & C) - B), len((B & C) - A),
        len(A & B & C),
    )

def make_layout_subsets(A, B, C):
    sA  = max(1, len(A)) ** 0.5
    sB  = max(1, len(B)) ** 0.5
    sC  = max(1, len(C)) ** 0.5
    ov  = min(sA, sB, sC) * 0.18
    return (
        max(0.05, sA - 2*ov), max(0.05, sB - 2*ov), ov,
        max(0.05, sC - 2*ov), ov, ov, ov * 0.4,
    )

def draw_venn(ax, A, B, C, set_labels, color, label_size=28):
    actual = subset_sizes(A, B, C)
    layout = make_layout_subsets(A, B, C)
    v = venn3(subsets=layout, set_labels=set_labels, ax=ax)

    for lbl, val in zip(v.subset_labels, actual):
        if lbl is not None:
            lbl.set_text(str(val))
            lbl.set_fontsize(label_size - 2)
            lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    for pid, alpha in [('100',0.50),('010',0.50),('001',0.50),
                        ('110',0.65),('101',0.65),('011',0.65),('111',0.85)]:
        patch = v.get_patch_by_id(pid)
        if patch:
            patch.set_facecolor(color); patch.set_alpha(alpha); patch.set_edgecolor('none')

    vc = venn3_circles(subsets=layout, ax=ax)
    for circle in (vc.patches if hasattr(vc, 'patches') else vc):
        circle.set_edgecolor('#222222'); circle.set_linewidth(5.0); circle.set_facecolor('none')

    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(label_size); lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold');   lbl.set_color('#111111')
    return v

# ── Load all data upfront ──────────────────────────────────────────────────────
compartment_sets = []
for label, excel_name, prefix, fmode in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    n_cols  = [c for c in df.columns if any(prefix+s in str(c) for s in N_SUFFIXES)]
    i_cols  = [c for c in df.columns if any(prefix+s in str(c) for s in I_SUFFIXES)]
    aw_cols = [c for c in df.columns if any(prefix+s in str(c) for s in AW_SUFFIXES)]
    pa_cols = [c for c in df.columns if any(prefix+s in str(c) for s in PA_SUFFIXES)]

    aw_up = get_sig_genes(df, n_cols, aw_cols, 'up')
    i_up  = get_sig_genes(df, n_cols, i_cols,  'up')
    pa_up = get_sig_genes(df, n_cols, pa_cols, 'up')
    aw_dn = get_sig_genes(df, n_cols, aw_cols, 'down')
    i_dn  = get_sig_genes(df, n_cols, i_cols,  'down')
    pa_dn = get_sig_genes(df, n_cols, pa_cols, 'down')

    compartment_sets.append(dict(
        label=label,
        aw_up=aw_up, i_up=i_up, pa_up=pa_up,
        aw_dn=aw_dn, i_dn=i_dn, pa_dn=pa_dn,
    ))
    print(f'{label}: AW↑={len(aw_up)} AW↓={len(aw_dn)}  (AW-M-3 excluded)')

# ── Individual pages (one per compartment) ────────────────────────────────────
with PdfPages(OUT_INDIVIDUAL) as pdf:
    for d in compartment_sets:
        label = d['label']
        fig, axes = plt.subplots(1, 2, figsize=(36, 16))
        fig.patch.set_facecolor('white')

        draw_venn(axes[0], d['aw_up'], d['i_up'], d['pa_up'],
                  set_labels=(f"AW Increased\n(n={len(d['aw_up'])})",
                               f"Intox. Increased\n(n={len(d['i_up'])})",
                               f"PA Increased\n(n={len(d['pa_up'])})"),
                  color=C_UP, label_size=28)
        axes[0].set_title(f'{label} — Proteins Increased',
                          fontsize=34, fontweight='bold', fontfamily='Arial', pad=20)

        draw_venn(axes[1], d['aw_dn'], d['i_dn'], d['pa_dn'],
                  set_labels=(f"AW Decreased\n(n={len(d['aw_dn'])})",
                               f"Intox. Decreased\n(n={len(d['i_dn'])})",
                               f"PA Decreased\n(n={len(d['pa_dn'])})"),
                  color=C_DOWN, label_size=28)
        axes[1].set_title(f'{label} — Proteins Decreased',
                          fontsize=34, fontweight='bold', fontfamily='Arial', pad=20)

        fig.suptitle(
            f'{label}  —  AW vs Naïve overlap with Intoxication & Protracted Abstinence\n'
            f'|FC| > {FC_THRESH}, corrected p > {CORR_THRESH}  |  AW-M-3 excluded',
            fontsize=26, fontfamily='Arial', color='#444444', style='italic', y=0.02
        )
        plt.tight_layout(pad=3.0, rect=[0, 0.06, 1, 1])
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)

print(f'Saved: {OUT_INDIVIDUAL}')

# ── Summary page — all 8 Venns on one page (4 rows × 2 cols) ─────────────────
fig, axes = plt.subplots(4, 2, figsize=(32, 48))
fig.patch.set_facecolor('white')

for row, d in enumerate(compartment_sets):
    label = d['label']

    # Left col = Increased
    draw_venn(axes[row, 0], d['aw_up'], d['i_up'], d['pa_up'],
              set_labels=(f"AW↑\n(n={len(d['aw_up'])})",
                           f"Intox.↑\n(n={len(d['i_up'])})",
                           f"PA↑\n(n={len(d['pa_up'])})"),
              color=C_UP, label_size=18)
    axes[row, 0].set_title(f'{label} — Increased',
                            fontsize=22, fontweight='bold', fontfamily='Arial', pad=14)

    # Right col = Decreased
    draw_venn(axes[row, 1], d['aw_dn'], d['i_dn'], d['pa_dn'],
              set_labels=(f"AW↓\n(n={len(d['aw_dn'])})",
                           f"Intox.↓\n(n={len(d['i_dn'])})",
                           f"PA↓\n(n={len(d['pa_dn'])})"),
              color=C_DOWN, label_size=18)
    axes[row, 1].set_title(f'{label} — Decreased',
                            fontsize=22, fontweight='bold', fontfamily='Arial', pad=14)

fig.suptitle(
    'Overlap of Significantly Changed Proteins Across Conditions\n'
    f'|FC| > {FC_THRESH}, corrected p > {CORR_THRESH}  |  AW-M-3 excluded',
    fontsize=28, fontweight='bold', fontfamily='Arial', y=0.995
)
plt.tight_layout(pad=3.0, rect=[0, 0, 1, 0.993])

with PdfPages(OUT_SUMMARY) as pdf:
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved: {OUT_SUMMARY}')
