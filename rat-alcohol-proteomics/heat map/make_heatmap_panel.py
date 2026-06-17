"""
Heatmap panel — all 4 compartments stacked vertically.
- 3 rows per compartment (Intoxication, Acute Withdrawal, Protracted Abstinence)
- Proteins as columns, sorted by AW fold change
- Width of each band is proportional to number of proteins
- Color = Log2 Fold Change vs Naive  (blue = decreased, red = increased)
- AW-M-3 excluded from Acute Withdrawal
- Chromatin: Keep+Review filter.  All others: unfiltered.

Produces TWO PDFs:
  Heatmap_Panel_Union.pdf   — proteins significant in >=1 condition
  Heatmap_Panel_AllSig.pdf  — proteins significant in ALL 3 conditions
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

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# Blue → white → red  (symmetric diverging)
CMAP = LinearSegmentedColormap.from_list(
    'bwr',
    ['#2166AC', '#4393C3', '#92C5DE', '#D1E5F0',
     '#FFFFFF',
     '#FDDBC7', '#F4A582', '#D6604D', '#B2182B'],
    N=512
)

# Conditions (AW-M-3 excluded)
CONDITIONS = [
    ('Intoxication',          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2']),
    ('Acute Withdrawal',      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']),
    ('Protracted Abstinence', ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2']),
]
NAIVE_SUFFIXES = ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']
COND_LABELS    = [c for c, _ in CONDITIONS]

# Display order: top → bottom
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

# ── Load everything once ───────────────────────────────────────────────────────
preloaded = []
for (disp, sheet, prefix, fmode) in SHEETS:
    raw = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in raw.columns:
        df = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = raw.copy().reset_index(drop=True)

    naive_cols = [c for c in df.columns if any(prefix+s in str(c) for s in NAIVE_SUFFIXES)]

    cond_cols_map = {}
    fc_map   = {}
    corr_map = {}
    for cname, csuf in CONDITIONS:
        cc = [c for c in df.columns if any(prefix+s in str(c) for s in csuf)]
        cond_cols_map[cname] = cc
        if cc:
            fc, corr = calc_stats(df, naive_cols, cc)
            fc_map[cname]   = fc
            corr_map[cname] = corr

    union_mask  = pd.Series(False, index=df.index)
    allsig_mask = pd.Series(True,  index=df.index)
    for cname, _ in CONDITIONS:
        if cname in fc_map:
            sig = (corr_map[cname].notna() & fc_map[cname].notna() &
                   (corr_map[cname] > CORR_THRESH) & (fc_map[cname].abs() > FC_THRESH))
            union_mask  = union_mask  | sig
            allsig_mask = allsig_mask & sig

    preloaded.append(dict(disp=disp, df=df, naive_cols=naive_cols,
                          cond_cols_map=cond_cols_map,
                          union_mask=union_mask, allsig_mask=allsig_mask))
    print(f'{disp}: union={union_mask.sum()}  all-sig={allsig_mask.sum()}')

def build_panel(preloaded, mask_key):
    """Return list of dicts with FC matrix per compartment."""
    out = []
    for entry in preloaded:
        disp = entry['disp']
        df   = entry['df']
        mask = entry[mask_key]
        df_s = df[mask].reset_index(drop=True)
        n = len(df_s)
        if n == 0:
            print(f'  {disp}: 0 proteins, skipping')
            continue
        naive_v = df_s[entry['naive_cols']].apply(pd.to_numeric, errors='coerce')
        rows = []
        for cname, _ in CONDITIONS:
            cc = entry['cond_cols_map'][cname]
            if cc:
                cond_v = df_s[cc].apply(pd.to_numeric, errors='coerce')
                rows.append((cond_v.mean(axis=1) - naive_v.mean(axis=1)).values)
            else:
                rows.append(np.zeros(n))
        mat = np.vstack(rows)                   # (3, n_proteins)
        order = np.argsort(mat[1])              # sort by AW fold change
        out.append(dict(label=disp, n=n, mat=mat[:, order]))
    return out

def make_figure(panel_data, version_label, pt_per_prot=0.008, min_band_w=0.0, row_h=0.42):
    N_COMP = len(panel_data)
    N_COND = len(CONDITIONS)

    # ── Layout constants (all in inches) ─────────────────────────────────────
    PT_PER_PROT = pt_per_prot   # inches of width per protein column
    MIN_BAND_W  = min_band_w    # minimum band width (so tiny sets stay visible)
    ROW_H       = row_h         # height per condition row
    TITLE_H     = 0.38          # height above each band for compartment label
    BAND_GAP    = 0.12          # vertical gap between bottom of one band and title of next
    LEFT_IN     = 2.20          # left margin for compartment label text
    RIGHT_IN    = 2.60          # right margin for condition label text
    TOP_PAD     = 0.50          # above first title
    BOT_PAD     = 0.80          # below last band (colorbar lives here)
    CBAR_H      = 0.14          # colorbar height

    def band_w(n): return max(n * PT_PER_PROT, MIN_BAND_W)

    max_hm_w = max(band_w(d['n']) for d in panel_data)

    BAND_H = N_COND * ROW_H

    FIG_W = LEFT_IN + max_hm_w + RIGHT_IN
    FIG_H = (TOP_PAD
             + N_COMP * (TITLE_H + BAND_H)
             + (N_COMP - 1) * BAND_GAP
             + BOT_PAD)

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor('white')

    vmax = 2.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im_ref = None

    for idx, d in enumerate(panel_data):
        label  = d['label']
        n_prot = d['n']
        mat    = d['mat']   # (3, n_prot)
        hm_w   = band_w(n_prot)

        # y position: count from top downward
        y_title_top = FIG_H - TOP_PAD - idx * (TITLE_H + BAND_H + BAND_GAP)
        y_band_top  = y_title_top - TITLE_H
        y_band_bot  = y_band_top  - BAND_H

        # Convert to figure fractions
        ax_left   = LEFT_IN / FIG_W
        ax_bottom = y_band_bot / FIG_H
        ax_w      = hm_w / FIG_W
        ax_h      = BAND_H / FIG_H

        ax = fig.add_axes([ax_left, ax_bottom, ax_w, ax_h])
        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                       interpolation='nearest')
        im_ref = im

        # White dividers between condition rows
        for r in range(1, N_COND):
            ax.axhline(r - 0.5, color='white', linewidth=1.8, zorder=3)

        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        # ── Compartment label — left of each band, vertically centred ─────────
        fig.text(
            (LEFT_IN - 0.12) / FIG_W,
            (y_title_top - TITLE_H / 2) / FIG_H,
            f'{label}\n(n = {n_prot})',
            ha='right', va='center',
            fontsize=11, fontweight='bold', fontfamily='Arial', color='#111111'
        )

        # ── Condition labels — right of this band ─────────────────────────────
        for r, cname in enumerate(COND_LABELS):
            # centre of row r (rows go top→bottom in imshow, so row 0 is at top)
            y_row_center = y_band_top - (r + 0.5) * ROW_H
            fig.text(
                (LEFT_IN + hm_w + 0.15) / FIG_W,
                y_row_center / FIG_H,
                cname,
                ha='left', va='center',
                fontsize=9.5, fontfamily='Arial', color='#111111'
            )

        # Thin border around heatmap
        for sp in ['top','bottom','left','right']:
            ax.spines[sp].set_visible(True)
            ax.spines[sp].set_linewidth(0.5)
            ax.spines[sp].set_color('#AAAAAA')

    # ── Colorbar — centred at the bottom ──────────────────────────────────────
    cbar_w_in  = 2.5
    cbar_l_in  = (FIG_W - cbar_w_in) / 2
    cbar_bot   = (BOT_PAD * 0.15) / FIG_H
    ax_cbar    = fig.add_axes([cbar_l_in / FIG_W, cbar_bot,
                                cbar_w_in / FIG_W, CBAR_H / FIG_H])
    cb = fig.colorbar(im_ref, cax=ax_cbar, orientation='horizontal')
    cb.set_label('')          # no label — avoids overlap with ticks
    cb.ax.tick_params(labelsize=8)
    cb.set_ticks([-2, -1, 0, 1, 2])
    # Add a small caption below the colorbar instead
    ax_cbar.text(0.5, -0.9, 'Log\u2082 Fold Change vs Na\u00efve',
                 transform=ax_cbar.transAxes,
                 ha='center', va='top', fontsize=8, fontfamily='Arial', color='#444444')

    # footnote removed

    return fig

# ── Build and save ─────────────────────────────────────────────────────────────
# pt_per_prot : inches per protein column
# min_band_w  : minimum band width in inches (keeps tiny sets visible)
# row_h       : height per condition row in inches
RUNS = [
    ('union_mask',  'Heatmap_Panel_Union.pdf',
     'significant in \u22651 condition (corrected p > 3.3, |FC| > 0.5)',
     dict(pt_per_prot=0.008, min_band_w=0.0,  row_h=0.42)),

    ('allsig_mask', 'Heatmap_Panel_AllSig.pdf',
     'significant in ALL 3 conditions (corrected p > 3.3, |FC| > 0.5 each)',
     dict(pt_per_prot=0.08,  min_band_w=2.5,  row_h=0.65)),
]

for mask_key, out_pdf, note, layout in RUNS:
    panel_data = build_panel(preloaded, mask_key)
    fig = make_figure(panel_data, note, **layout)
    pdf = PdfPages(out_pdf)
    pdf.savefig(fig, dpi=200, bbox_inches='tight', facecolor='white')
    pdf.close()
    plt.close(fig)
    print(f'Saved: {out_pdf}')
