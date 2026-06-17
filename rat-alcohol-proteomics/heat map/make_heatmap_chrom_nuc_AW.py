"""
Heatmap: Chromatin vs Soluble Nuclear during Acute Withdrawal
- Only 2 compartments: Chromatin and Soluble Nuclear
- Proteins: union of proteins significant in AW in EITHER compartment
  (corrected p > 3.3, |FC| > 0.5 in AW vs Naïve)
- Color = Log2 Fold Change vs Naïve (same as Heatmap_Panel_Union.pdf)
- 3 condition rows: Intoxication, Acute Withdrawal, Protracted Abstinence
- Sorted by AW fold change (low → high)
- AW-M-3 excluded. Chromatin: Keep+Review filter.

Output: Heatmap_Chrom_vs_Nuclear_AW.pdf
"""

import pandas as pd
import numpy as np
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'Heatmap_Chrom_vs_Nuclear_AW.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

CMAP = LinearSegmentedColormap.from_list(
    'bwr',
    ['#2166AC', '#4393C3', '#92C5DE', '#D1E5F0',
     '#FFFFFF',
     '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
    N=512
)

COMPARTMENT_COLORS = {
    'Chromatin':       '#F7F3E3',
    'Soluble Nuclear': '#E3F2FA',
}

CONDITIONS = [
    ('Acute Withdrawal', ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']),
]
NAIVE_SUFFIXES = ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']
AW_SUFFIXES    = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']
COND_LABELS    = [c for c, _ in CONDITIONS]

SHEETS = [
    ('Chromatin',       'Chromatin',       'Chrom_', 'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'Nuc_',   'all'),
]

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

# ── Load both compartments ─────────────────────────────────────────────────────
panel_data = []

for (disp, sheet, prefix, fmode) in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    naive_cols = [c for c in df.columns if any(prefix+s in str(c) for s in NAIVE_SUFFIXES)]
    aw_cols    = [c for c in df.columns if any(prefix+s in str(c) for s in AW_SUFFIXES)]

    # Significant in AW only
    aw_fc, aw_corr = calc_stats(df, naive_cols, aw_cols)
    aw_sig = (aw_corr.notna() & aw_fc.notna() &
              (aw_corr > CORR_THRESH) & (aw_fc.abs() > FC_THRESH))

    df_sig = df[aw_sig].reset_index(drop=True)
    n = len(df_sig)
    print(f'{disp}: {n} proteins significant in AW')

    if n == 0:
        continue

    naive_v = df_sig[naive_cols].apply(pd.to_numeric, errors='coerce')

    # FC matrix: (3, n_proteins)
    rows = []
    for cname, csuf in CONDITIONS:
        cc = [c for c in df_sig.columns if any(prefix+s in str(c) for s in csuf)]
        if cc:
            cond_v = df_sig[cc].apply(pd.to_numeric, errors='coerce')
            rows.append((cond_v.mean(axis=1) - naive_v.mean(axis=1)).values)
        else:
            rows.append(np.zeros(n))
    mat = np.vstack(rows)

    # Sort by AW fold change (only 1 row so index 0)
    order = np.argsort(mat[0])
    mat   = mat[:, order]

    panel_data.append(dict(label=disp, n=n, mat=mat))

# ── Figure — same layout as Heatmap_Panel_Union ────────────────────────────────
N_COMP = len(panel_data)
N_COND = len(CONDITIONS)

FIG_W    = 18
ROW_H    = 0.55
BAND_GAP = 0.12
TITLE_H  = 0.38
TOP_PAD  = 0.45
BOT_PAD  = 0.80
LEFT_PAD = 2.20
RIGHT_PAD= 2.80

BAND_H = N_COND * ROW_H
fig_h  = TOP_PAD + N_COMP * (TITLE_H + BAND_H) + (N_COMP - 1) * BAND_GAP + BOT_PAD

fig = plt.figure(figsize=(FIG_W, fig_h))
fig.patch.set_facecolor('white')

vmax = 2.0
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im_ref = None

def ffy(y): return y / fig_h
def ffx(x): return x / FIG_W

for idx, d in enumerate(panel_data):
    label  = d['label']
    mat    = d['mat']
    n_prot = d['n']
    hm_w   = n_prot * 0.008

    band_top = fig_h - TOP_PAD - idx * (TITLE_H + BAND_H + BAND_GAP) - TITLE_H
    band_bot = band_top - BAND_H

    ax = fig.add_axes([ffx(LEFT_PAD), ffy(band_bot), ffx(hm_w), ffy(BAND_H)])
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    im_ref = im

    for r in range(1, N_COND):
        ax.axhline(r - 0.5, color='white', linewidth=1.5, zorder=3)

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

    # Condition labels on right
    for r, cname in enumerate(COND_LABELS):
        y_frac = ffy(band_bot + BAND_H - (r + 0.5) * ROW_H)
        fig.text(ffx(LEFT_PAD + hm_w + 0.15), y_frac, cname,
                 ha='left', va='center', fontsize=10,
                 fontfamily='Arial', color='#111111')

    # Compartment title + n=
    fig.text(ffx(LEFT_PAD - 0.1), ffy(band_top + TITLE_H * 0.55),
             label, ha='left', va='center',
             fontsize=13, fontweight='bold', fontfamily='Arial', color='#111111')
    fig.text(ffx(LEFT_PAD - 0.1), ffy(band_top + TITLE_H * 0.10),
             f'n = {n_prot} proteins', ha='left', va='center',
             fontsize=8.5, fontfamily='Arial', color='#555555', style='italic')

# ── Colorbar ───────────────────────────────────────────────────────────────────
ax_cbar = fig.add_axes([ffx((FIG_W - 2.5) / 2), ffy(0.22), ffx(2.5), ffy(0.18)])
cb = fig.colorbar(im_ref, cax=ax_cbar, orientation='horizontal')
cb.set_label('')
cb.ax.tick_params(labelsize=8)
cb.set_ticks([-2, -1, 0, 1, 2])
ax_cbar.text(0.5, -0.9, 'Log\u2082 Fold Change vs Na\u00efve',
             transform=ax_cbar.transAxes,
             ha='center', va='top', fontsize=8.5,
             fontfamily='Arial', color='#444444')

pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
pdf.close()
plt.close(fig)
print(f'Saved: {OUT}')
