import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Venn_AcuteWithdrawal.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# Bolder colours matching the VP panels
C_UP   = '#F28FAD'   # bold pink  (increased)
C_DOWN = '#7BAFD4'   # bold blue  (decreased)

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_N-',  'Memb_I-',  'Memb_AW-',  'Memb_PA-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-',  'Cyto_I-',  'Cyto_AW-',  'Cyto_PA-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_N-', 'Chrom_I-', 'Chrom_AW-', 'Chrom_PA-', 'keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_N-',   'Nuc_I-',   'Nuc_AW-',   'Nuc_PA-',   'keep_review'),
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

def get_sig_genes(df, naive_cols, cond_cols, direction):
    fc, corr = calc_stats(df, naive_cols, cond_cols)
    mask = fc.notna() & corr.notna() & (corr > CORR_THRESH)
    if direction == 'up':
        sig = mask & (fc > FC_THRESH)
    else:
        sig = mask & (fc < -FC_THRESH)
    return set(df.loc[sig, 'Gene symbol'].astype(str).tolist())

def subset_sizes(A, B, C):
    """Return tuple (Abc, aBc, ABc, abC, AbC, aBC, ABC) — actual counts."""
    return (
        len(A - B - C),
        len(B - A - C),
        len((A & B) - C),
        len(C - A - B),
        len((A & C) - B),
        len((B & C) - A),
        len(A & B & C),
    )

def make_layout_subsets(A, B, C):
    """Return fake subsets that give circle sizes mildly proportional to set sizes.
    Uses n^0.5 scaling so the largest set is visibly but not dramatically bigger.
    All 7 regions are kept > 0 so every label is accessible."""
    sA = max(1, len(A)) ** 0.5
    sB = max(1, len(B)) ** 0.5
    sC = max(1, len(C)) ** 0.5
    ov = min(sA, sB, sC) * 0.18   # small equal overlap for all intersection regions
    return (
        max(0.05, sA - 2 * ov),    # Abc  (A only)
        max(0.05, sB - 2 * ov),    # aBc  (B only)
        ov,                         # ABc  (A∩B only)
        max(0.05, sC - 2 * ov),    # abC  (C only)
        ov,                         # AbC  (A∩C only)
        ov,                         # aBC  (B∩C only)
        ov * 0.4,                   # ABC  (all three)
    )

def draw_venn(ax, A, B, C, set_labels, color, label_size=28):
    """Draw a 3-way Venn: largest set slightly bigger; all region counts always shown."""
    actual = subset_sizes(A, B, C)
    layout = make_layout_subsets(A, B, C)

    v = venn3(subsets=layout, set_labels=set_labels, ax=ax)

    # Overwrite every region label with the actual count (always show, even 0)
    for lbl, val in zip(v.subset_labels, actual):
        if lbl is not None:
            lbl.set_text(str(val))
            lbl.set_fontsize(label_size - 2)
            lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    # Colour the patches
    for patch_id, alpha in [
        ('100', 0.50), ('010', 0.50), ('001', 0.50),
        ('110', 0.65), ('101', 0.65), ('011', 0.65),
        ('111', 0.85),
    ]:
        patch = v.get_patch_by_id(patch_id)
        if patch:
            patch.set_facecolor(color)
            patch.set_alpha(alpha)
            patch.set_edgecolor('none')

    # Draw thick circle outlines on top (same layout so circles match)
    vc = venn3_circles(subsets=layout, ax=ax)
    for circle in (vc.patches if hasattr(vc, 'patches') else vc):
        circle.set_edgecolor('#222222')
        circle.set_linewidth(5.0)
        circle.set_facecolor('none')

    # Style set labels
    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(label_size)
            lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    return v

with PdfPages(OUT) as pdf:
    for label, excel_name, n_pat, i_pat, aw_pat, pa_pat, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
        if fmode == 'keep_review':
            df = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        naive_cols = [c for c in df.columns if n_pat  in str(c)]
        i_cols     = [c for c in df.columns if i_pat  in str(c)]
        aw_cols    = [c for c in df.columns if aw_pat in str(c)]
        pa_cols    = [c for c in df.columns if pa_pat in str(c)]

        # Sets for Venn 1 (increased in AW, Intox, and PA)
        aw_up = get_sig_genes(df, naive_cols, aw_cols, 'up')
        i_up  = get_sig_genes(df, naive_cols, i_cols,  'up')
        pa_up = get_sig_genes(df, naive_cols, pa_cols, 'up')

        # Sets for Venn 2 (decreased in AW, Intox, and PA)
        aw_dn = get_sig_genes(df, naive_cols, aw_cols, 'down')
        i_dn  = get_sig_genes(df, naive_cols, i_cols,  'down')
        pa_dn = get_sig_genes(df, naive_cols, pa_cols, 'down')

        fig, axes = plt.subplots(1, 2, figsize=(36, 16))
        fig.patch.set_facecolor('white')

        # ── Venn 1: increased in AW, Intox, and PA ──
        draw_venn(
            axes[0], aw_up, i_up, pa_up,
            set_labels=(
                f'AW Increased\n(n={len(aw_up)})',
                f'Intox. Increased\n(n={len(i_up)})',
                f'PA Increased\n(n={len(pa_up)})',
            ),
            color=C_UP, label_size=28
        )
        axes[0].set_title(f'{label} — Proteins Increased in AW',
                          fontsize=34, fontweight='bold', fontfamily='Arial', pad=20)

        # ── Venn 2: decreased in AW, Intox, and PA ──
        draw_venn(
            axes[1], aw_dn, i_dn, pa_dn,
            set_labels=(
                f'AW Decreased\n(n={len(aw_dn)})',
                f'Intox. Decreased\n(n={len(i_dn)})',
                f'PA Decreased\n(n={len(pa_dn)})',
            ),
            color=C_DOWN, label_size=28
        )
        axes[1].set_title(f'{label} — Proteins Decreased in AW',
                          fontsize=34, fontweight='bold', fontfamily='Arial', pad=20)

        fig.suptitle(
            f'{label}  —  Acute Withdrawal vs Naive\n'
            f'Overlap with Intoxication and Protracted Abstinence (|FC| > {FC_THRESH}, corrected p > {CORR_THRESH})',
            fontsize=30, fontfamily='Arial', color='#444444', style='italic', y=0.02
        )

        plt.tight_layout(pad=3.0, rect=[0, 0.06, 1, 1])
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'{label}: AW↑={len(aw_up)} AW↓={len(aw_dn)}')

print(f'Saved: {OUT}')
