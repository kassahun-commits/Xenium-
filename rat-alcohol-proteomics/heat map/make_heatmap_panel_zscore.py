"""
Heatmap panel — SAME visual style as Heatmap_Panel_Union.pdf but:
- Rows = individual replicates (not condition averages)
  Naïve x5, Intoxication x5, Acute Withdrawal x4 (AW-M-3 excluded), PA x5 = 19 rows
- Color = Z-score of LFQ intensity per protein across all 19 replicates
- Proteins sorted by AW mean fold change (low → high, same as union panel)
- Thin white dividers between condition groups
- Condition label spans its group of replicates on the right
- Chromatin: Keep+Review filter. All others: unfiltered.

Output: Heatmap_Panel_Zscore_Union.pdf
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
OUT  = 'Heatmap_Panel_Zscore_Union.pdf'

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
    'Membrane':        '#EDE8F5',
    'Cytosol':         '#E4F5E8',
    'Chromatin':       '#F7F3E3',
    'Soluble Nuclear': '#E3F2FA',
}

# Individual replicates grouped by condition — AW-M-3 excluded
GROUPS = [
    ('Naïve',                 ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']),
    ('Intoxication',          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2']),
    ('Acute Withdrawal',      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']),
    ('Protracted Abstinence', ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2']),
]
NAIVE_SUFFIXES = ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']
AW_SUFFIXES    = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_',  'all'),
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

# ── Load data ──────────────────────────────────────────────────────────────────
panel_data = []

for (disp, sheet, prefix, fmode) in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    naive_cols = [c for c in df.columns if any(prefix+s in str(c) for s in NAIVE_SUFFIXES)]
    aw_cols    = [c for c in df.columns if any(prefix+s in str(c) for s in AW_SUFFIXES)]

    # Build all replicate column lists per group
    group_rep_cols = []
    all_rep_cols   = []
    for gname, gsuf in GROUPS:
        gc = [c for c in df.columns if any(prefix+s in str(c) for s in gsuf)]
        group_rep_cols.append((gname, gc))
        all_rep_cols.extend(gc)

    # Union sig mask
    union_mask = pd.Series(False, index=df.index)
    for gname, gsuf in GROUPS:
        if gname == 'Naïve': continue
        gc = [c for c in df.columns if any(prefix+s in str(c) for s in gsuf)]
        if gc:
            fc, corr = calc_stats(df, naive_cols, gc)
            sig = (corr.notna() & fc.notna() &
                   (corr > CORR_THRESH) & (fc.abs() > FC_THRESH))
            union_mask = union_mask | sig

    df_s = df[union_mask].reset_index(drop=True)
    n = len(df_s)
    print(f'{disp}: {n} proteins')
    if n == 0:
        continue

    # Sort by AW mean fold change
    naive_v = df_s[naive_cols].apply(pd.to_numeric, errors='coerce')
    aw_v    = df_s[aw_cols].apply(pd.to_numeric, errors='coerce')
    aw_fc   = (aw_v.mean(axis=1) - naive_v.mean(axis=1)).values
    order   = np.argsort(aw_fc)
    df_s    = df_s.iloc[order].reset_index(drop=True)

    # Z-score each protein across all 19 replicates
    mat_raw = df_s[all_rep_cols].apply(pd.to_numeric, errors='coerce').values.astype(float)
    # mat_raw shape: (n_proteins, n_reps)
    row_mean = np.nanmean(mat_raw, axis=1, keepdims=True)
    row_std  = np.nanstd(mat_raw,  axis=1, keepdims=True)
    row_std[row_std == 0] = 1
    mat_z = (mat_raw - row_mean) / row_std  # (n_proteins, n_reps)

    # Transpose so rows=reps, cols=proteins  → imshow(mat) gives correct orientation
    mat_plot = mat_z.T   # (n_reps, n_proteins)

    panel_data.append({
        'label':          disp,
        'n':              n,
        'mat':            mat_plot,        # (n_reps, n_proteins)
        'group_rep_cols': group_rep_cols,  # [(gname, [cols]), ...]
    })

# ── Figure layout — same as Heatmap_Panel_Union ────────────────────────────────
N_COMP   = len(panel_data)
N_REPS   = sum(len(gc) for _, gc in panel_data[0]['group_rep_cols'])   # 19
N_GROUPS = len(GROUPS)

FIG_W    = 18
ROW_H    = 0.20    # inches per replicate row (thinner than condition rows)
BAND_GAP = 0.12
TITLE_H  = 0.38
TOP_PAD  = 0.45
BOT_PAD  = 0.80
LEFT_PAD = 2.20
RIGHT_PAD= 2.80

BAND_H = N_REPS * ROW_H

fig_h = TOP_PAD + N_COMP * (TITLE_H + BAND_H) + (N_COMP - 1) * BAND_GAP + BOT_PAD

fig = plt.figure(figsize=(FIG_W, fig_h))
fig.patch.set_facecolor('white')

vmax = 2.5
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im_ref = None

def ffy(y): return y / fig_h
def ffx(x): return x / FIG_W

for idx, d in enumerate(panel_data):
    label          = d['label']
    mat            = d['mat']            # (n_reps, n_proteins)
    n_prot         = d['n']
    group_rep_cols = d['group_rep_cols']
    hm_w           = n_prot * 0.008

    band_top = fig_h - TOP_PAD - idx * (TITLE_H + BAND_H + BAND_GAP) - TITLE_H
    band_bot = band_top - BAND_H

    ax = fig.add_axes([ffx(LEFT_PAD), ffy(band_bot), ffx(hm_w), ffy(BAND_H)])
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    im_ref = im

    # White dividers between condition GROUPS (thicker) and between reps (thinner)
    rep_cursor = 0
    for gname, gc in group_rep_cols:
        n_reps = len(gc)
        # Thick white line between groups
        if rep_cursor > 0:
            ax.axhline(rep_cursor - 0.5, color='white', linewidth=2.5, zorder=3)
        # Thin lines between individual reps within group
        for r in range(1, n_reps):
            ax.axhline(rep_cursor + r - 0.5, color='white', linewidth=0.5, zorder=2)
        rep_cursor += n_reps

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

    # Condition group labels on right — centred on their group of rows
    rep_cursor = 0
    for gname, gc in group_rep_cols:
        n_reps = len(gc)
        # Centre of this group in figure fraction
        grp_center_y = ffy(band_bot + BAND_H - (rep_cursor + n_reps / 2) * ROW_H)
        fig.text(ffx(LEFT_PAD + hm_w + 0.15), grp_center_y, gname,
                 ha='left', va='center', fontsize=10,
                 fontfamily='Arial', color='#111111')
        rep_cursor += n_reps

    # Compartment title + n=
    fig.text(ffx(LEFT_PAD - 0.1), ffy(band_top + TITLE_H * 0.55),
             label, ha='left', va='center',
             fontsize=13, fontweight='bold', fontfamily='Arial', color='#111111')
    fig.text(ffx(LEFT_PAD - 0.1), ffy(band_top + TITLE_H * 0.10),
             f'n = {n_prot} proteins', ha='left', va='center',
             fontsize=8.5, fontfamily='Arial', color='#555555', style='italic')

# ── Colorbar ───────────────────────────────────────────────────────────────────
ax_cbar = fig.add_axes([ffx((FIG_W - 2.5) / 2), ffy(0.18), ffx(2.5), ffy(0.14)])
cb = fig.colorbar(im_ref, cax=ax_cbar, orientation='horizontal')
cb.set_label('')
cb.ax.tick_params(labelsize=8)
cb.set_ticks([-2, -1, 0, 1, 2])
ax_cbar.text(0.5, -0.9, 'Z-score (LFQ intensity)',
             transform=ax_cbar.transAxes,
             ha='center', va='top', fontsize=8.5,
             fontfamily='Arial', color='#444444')

pdf = PdfPages(OUT)
pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
pdf.close()
plt.close(fig)
print(f'Saved: {OUT}')
