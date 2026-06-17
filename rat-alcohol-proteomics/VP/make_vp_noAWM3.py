"""
Volcano plots for all 4 compartments — AW vs Naïve — with AW-M-3 EXCLUDED.
AW-M-3 is a global outlier replicate (median 2-3x lower than all other AW replicates
across every compartment), consistent with a failed sample run rather than biology.

Saves: VP_AW_noAWM3_AllCompartments.pdf   (new version, original untouched)
Also saves: AW_M3_outlier_check.pdf        (diagnostic plot showing the outlier)
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

FILE   = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_VP = 'VP_AW_noAWM3_AllCompartments.pdf'
OUT_DX = 'AW_M3_outlier_check.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_OTHER = '#C8C8C8'
C_NS    = '#EBEBEB'

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'Memb_N-',  ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'], 'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'Cyto_N-',  ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'], 'all'),
    ('Chromatin',       'Chromatin',       'Chrom_', 'Chrom_N-', ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'], 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'Nuc_N-',   ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'], 'all'),
]

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    pv    = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        pv.append(stats.ttest_ind(n, c, equal_var=False)[1]
                  if len(n) >= 2 and len(c) >= 2 else np.nan)
    pv    = pd.Series(pv, index=df.index)
    valid = pv.notna()
    ranks = pv[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pv[valid] * valid.sum() / ranks))
    return fc, corr

def draw_volcano(ax, fc, corrected, highlight, title):
    mask = fc.notna() & corrected.notna()
    f, cr = fc[mask].values, corrected[mask].values
    up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
    down = (cr > CORR_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down
    hi, lo, c_hi = (up, down, C_UP) if highlight == 'up' else (down, up, C_DOWN)

    ax.scatter(f[ns], cr[ns], s=12, color=C_NS,    alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(f[lo], cr[lo], s=22, color=C_OTHER, alpha=0.50, linewidths=0, rasterized=True)
    ax.scatter(f[hi], cr[hi], s=40, color=c_hi,    alpha=0.90, linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.8, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.8, alpha=0.4)

    ax.set_xlabel('Fold Change',                 fontsize=18, labelpad=8)
    ax.set_ylabel('Corrected p-value\n(−log₂)',  fontsize=18, labelpad=8)
    ax.set_title(title, fontsize=19, fontweight='bold', pad=10)
    ax.tick_params(labelsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_hi  = int(hi.sum())
    n_tot = int(mask.sum())
    direction = '↑' if highlight == 'up' else '↓'
    ax.text(0.98, 0.98, f'{direction}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=17, color='#333333', fontweight='bold')

    dir_label = 'Increased' if highlight == 'up' else 'Decreased'
    ax.legend(handles=[
        mpatches.Patch(color=c_hi,    label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_OTHER, label='Other sig.'),
        mpatches.Patch(color=C_NS,    label='NS'),
    ], fontsize=12, loc='upper left', framealpha=0.88,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.8)


# ── 1. Outlier diagnostic plot ─────────────────────────────────────────────────
print('Making outlier diagnostic plot...')
with PdfPages(OUT_DX) as pdf:
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    fig.patch.set_facecolor('white')
    fig.suptitle('AW Replicate Distributions — AW-M-3 Highlighted as Outlier',
                 fontsize=16, fontweight='bold', y=0.98)

    col = 0
    for disp_label, sheet, prefix, n_pat, aw_suf, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=sheet)
        if fmode == 'keep_review' and 'Filter' in raw.columns:
            df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        aw_cols    = [c for c in df.columns if any(prefix+s in str(c) for s in aw_suf)]
        naive_cols = [c for c in df.columns if n_pat in str(c)]

        # Row 1: boxplot of all replicates
        ax1 = axes[0, col]
        data_to_plot = []
        labels       = []
        colors       = []
        for c in naive_cols + aw_cols:
            v = pd.to_numeric(df[c], errors='coerce').dropna().values
            data_to_plot.append(v)
            lbl = c.replace(prefix, '').replace('Chrom_','').replace('Memb_','').replace('Cyto_','').replace('Nuc_','')
            labels.append(lbl)
            colors.append('#FFAAAA' if 'AW-M-3' in c else
                          ('#FFD580' if 'AW-' in c else '#AACCFF'))

        bp = ax1.boxplot(data_to_plot, patch_artist=True, showfliers=False,
                         medianprops=dict(color='black', linewidth=2))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        ax1.set_xticks(range(1, len(labels)+1))
        ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax1.set_title(f'{disp_label}\nReplicate distributions', fontsize=11, fontweight='bold')
        ax1.set_ylabel('LFQ intensity', fontsize=9)
        ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)

        # Row 2: median of each AW replicate as bar chart
        ax2 = axes[1, col]
        aw_medians = []
        aw_labels  = []
        bar_colors = []
        for c in aw_cols:
            v = pd.to_numeric(df[c], errors='coerce').dropna()
            aw_medians.append(v.median())
            aw_labels.append(c.replace(prefix,''))
            bar_colors.append('#E8305A' if 'AW-M-3' in c else '#2B7FD4')

        bars = ax2.bar(range(len(aw_medians)), aw_medians,
                       color=bar_colors, edgecolor='white', linewidth=0.5)
        ax2.set_xticks(range(len(aw_labels)))
        ax2.set_xticklabels(aw_labels, rotation=45, ha='right', fontsize=8)
        ax2.set_title('AW replicate medians\n(red = AW-M-3 outlier)', fontsize=11)
        ax2.set_ylabel('Median LFQ', fontsize=9)
        ax2.axhline(np.mean(aw_medians[:-1]), color='#555555',
                    linestyle='--', linewidth=1.5, label='Mean excl. AW-M-3')
        ax2.legend(fontsize=8)
        ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)

        col += 1

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close(fig)

print(f'Saved: {OUT_DX}')


# ── 2. Volcano plots WITHOUT AW-M-3 ───────────────────────────────────────────
print('Making volcano plots without AW-M-3...')
AW_SUBS = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']   # AW-M-3 excluded

with PdfPages(OUT_VP) as pdf:
    for disp_label, sheet, prefix, n_pat, _, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=sheet)
        if fmode == 'keep_review' and 'Filter' in raw.columns:
            df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
        else:
            df = raw.copy().reset_index(drop=True)

        naive_cols = [c for c in df.columns if n_pat in str(c)]
        aw_cols_4  = [c for c in df.columns if any(prefix+s in str(c) for s in AW_SUBS)]

        fc, corrected = calc_stats(df, naive_cols, aw_cols_4)

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.patch.set_facecolor('white')

        draw_volcano(axes[0], fc, corrected, 'up',
                     f'{disp_label}\nAW vs Naïve — Increased  (excl. AW-M-3)')
        draw_volcano(axes[1], fc, corrected, 'down',
                     f'{disp_label}\nAW vs Naïve — Decreased  (excl. AW-M-3)')

        fig.text(0.5, 0.01,
                 f'Thresholds: |FC| > {FC_THRESH}  |  Corrected p > {CORR_THRESH}'
                 f'  |  AW replicates: AW-F-1, AW-F-2, AW-M-1, AW-M-2  (AW-M-3 excluded)',
                 ha='center', fontsize=10, color='#666666', style='italic')

        plt.tight_layout(pad=2.0, rect=[0, 0.04, 1, 1])
        pdf.savefig(fig, dpi=180, bbox_inches='tight')
        plt.close(fig)

        n_up   = int(((corrected > CORR_THRESH) & (fc >  FC_THRESH)).sum())
        n_down = int(((corrected > CORR_THRESH) & (fc < -FC_THRESH)).sum())
        print(f'  {disp_label}: ↑{n_up}  ↓{n_down}')

print(f'Saved: {OUT_VP}')
