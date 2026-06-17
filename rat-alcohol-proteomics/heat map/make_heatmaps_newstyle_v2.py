"""
Heatmaps — publication style v2
================================
Remakes three heatmaps in a clean publication aesthetic.
All bands share the same fixed width so the figure fits on a page.
Overrides the _v2 files only — originals untouched.

Outputs:
  Heatmap_Union_PubQuality_v2.pdf/.png
  Heatmap_Panel_Union_v2.pdf/.png
  Heatmap_Chrom_vs_Nuclear_AW_v2.pdf/.png
"""

import os, shutil
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE = os.path.join(SCRIPT_DIR,
       '../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')

# ── Colormap ──────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0','#2B7FD4','#6AAEE0','#AECFE8','#DDEEF8',
     '#FFFFFF',
     '#FDD8E7','#F5A0BC','#EE5F8B','#E8305A','#C01E42'],
    N=512,
)
VMAX = 3.0
norm = TwoSlopeNorm(vmin=-VMAX, vcenter=0, vmax=VMAX)

# ── Page layout (proportional band widths, figure width set by widest band) ───
INCHES_PER_PROT = 0.004   # width per protein column — Chromatin (~1200) ≈ 4.9 in
LEFT_IN         = 1.60    # left margin for compartment labels
RIGHT_IN        = 2.00    # right margin for condition labels

ROW_H       = 0.26    # inches per condition row (shorter = tighter)
N_COND      = 3
BAND_H      = N_COND * ROW_H            # 0.78 inches per compartment
SECTION_GAP = 0.18    # gap between compartments
TOP_PAD     = 0.45
BOT_PAD     = 0.95    # room for colorbar + footnote
CBAR_W_IN   = 2.6
CBAR_H_IN   = 0.12
CBAR_BOT    = 0.38

CORR_THRESH = 3.3
FC_THRESH   = 0.5

COND_LABELS = ['Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']

COMP_BG = {
    'Membrane':        '#F7F3FA',
    'Cytosol':         '#F2F9F4',
    'Chromatin':       '#FBF9F0',
    'Soluble Nuclear': '#F0F6FB',
}

NAIVE_SUFFIXES = ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2']
INTOX_SUFFIXES = ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2']
AW_SUFFIXES    = ['AW-F-1','AW-F-2','AW-M-1','AW-M-2']
PA_SUFFIXES    = ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2']

# ── Stat helpers ──────────────────────────────────────────────────────────────
def get_cols(df, prefix, suffixes):
    return [c for c in df.columns if any(prefix + s in str(c) for s in suffixes)]

def mean_fc(df, naive_cols, cond_cols):
    n = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    c = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    return c.mean(axis=1) - n.mean(axis=1)

def corrected_p(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    pvs = []
    for i in range(len(df)):
        nv = naive.iloc[i].dropna().values
        cv = cond.iloc[i].dropna().values
        pvs.append(stats.ttest_ind(nv, cv, equal_var=False)[1]
                   if len(nv) >= 2 and len(cv) >= 2 else np.nan)
    pvs   = pd.Series(pvs, index=df.index)
    valid = pvs.notna()
    ranks = pvs[valid].rank(ascending=False)
    corr  = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, pvs[valid] * valid.sum() / ranks))
    return corr

# ── Core draw function ────────────────────────────────────────────────────────
def draw(panel_data, title, footnote, out_pdf, out_png, n_cond=3, cond_labels=None):
    """
    panel_data : list of dicts — label, n, mat (n_cond × n_prot)
    Each band width is proportional to its protein count.
    Figure width is set by the widest band (Chromatin).
    """
    if cond_labels is None:
        cond_labels = COND_LABELS[:n_cond]
    # Sort least → most proteins (fewest at top)
    panel_data = sorted(panel_data, key=lambda d: d['n'])
    band_h  = n_cond * ROW_H
    max_hm_w = max(d['n'] * INCHES_PER_PROT for d in panel_data)
    FIG_W   = LEFT_IN + max_hm_w + RIGHT_IN

    FIG_H = (TOP_PAD
             + len(panel_data) * band_h
             + (len(panel_data) - 1) * SECTION_GAP
             + BOT_PAD)

    fig      = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    im_ref   = None
    y_cursor = FIG_H - TOP_PAD

    for d in panel_data:
        label  = d['label']
        n_prot = d['n']
        mat    = d['mat']
        hm_w   = n_prot * INCHES_PER_PROT   # proportional width for this band

        y_bot = y_cursor - band_h

        ax = fig.add_axes([LEFT_IN / FIG_W, y_bot / FIG_H,
                           hm_w / FIG_W,   band_h / FIG_H],
                          facecolor=COMP_BG.get(label, '#F8F8F8'))

        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                       interpolation='nearest')
        if im_ref is None:
            im_ref = im

        # Row dividers
        for r in range(1, n_cond):
            ax.axhline(r - 0.5, color='white', linewidth=1.0, zorder=3)

        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

        # Condition labels — right side
        for r, cname in enumerate(cond_labels):
            y_frac = (y_bot + band_h - (r + 0.5) * ROW_H) / FIG_H
            fig.text((LEFT_IN + hm_w + 0.10) / FIG_W, y_frac,
                     cname, ha='left', va='center',
                     fontsize=7.5, fontfamily='Arial', color='#222222')

        # Compartment label — left side
        mid_y = (y_bot + band_h / 2) / FIG_H
        fig.text((LEFT_IN - 0.10) / FIG_W, mid_y,
                 label, ha='right', va='center',
                 fontsize=9, fontweight='bold',
                 fontfamily='Arial', color='#111111')

        # n= in smaller italic below
        fig.text((LEFT_IN - 0.10) / FIG_W, mid_y - 0.048,
                 f'n = {n_prot}', ha='right', va='top',
                 fontsize=6.5, fontfamily='Arial',
                 color='#666666', style='italic')

        y_cursor = y_bot - SECTION_GAP

    # Colorbar
    cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
    ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H,
                           CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
    cb = fig.colorbar(im_ref, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
    cb.ax.tick_params(labelsize=6.5)
    ax_cb.set_xlabel('Log2 Fold Change vs Naive',
                     fontsize=7.5, fontfamily='Arial',
                     color='#444444', labelpad=3)

    # Footnote
    fig.text(0.5, 0.01, footnote,
             ha='center', va='bottom', fontsize=5,
             fontfamily='Arial', color='#999999')

    # Title
    fig.text(0.5, (FIG_H - TOP_PAD * 0.32) / FIG_H,
             title, ha='center', va='top',
             fontsize=8.5, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_png}')


# ══════════════════════════════════════════════════════════════════════════════
# HEATMAP 1 — Union PubQuality (pre-computed FC columns)
# ══════════════════════════════════════════════════════════════════════════════
print('\n=== Heatmap 1: Union PubQuality v2 ===')
COMPARTMENTS_4 = [
    ('Membrane',        'Membrane',        'all'),
    ('Cytosol',         'Cytosol',         'all'),
    ('Chromatin',       'Chromatin',       'keep_review'),
    ('Soluble Nuclear', 'Soluble nuclear', 'all'),
]

panel_data_pub = []
for disp, sheet, fmode in COMPARTMENTS_4:
    df = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in df.columns:
        df = df[df['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    union_mask = pd.Series(False, index=df.index)
    for col in ['Corrected', 'Corrected.1', 'Corrected.2']:
        if col in df.columns:
            union_mask |= pd.to_numeric(df[col], errors='coerce').ge(CORR_THRESH).fillna(False)

    df_u = df[union_mask].reset_index(drop=True)
    n    = len(df_u)
    print(f'  {disp}: {n}')
    if n == 0: continue

    rows = []
    for fc_col in ['Fold change', 'Fold change.1', 'Fold change.2']:
        rows.append(pd.to_numeric(df_u.get(fc_col, pd.Series(np.zeros(n))),
                                  errors='coerce').fillna(0).values)
    mat   = np.vstack(rows)
    order = np.argsort(mat[1])[::-1]
    panel_data_pub.append(dict(label=disp, n=n, mat=mat[:, order]))

draw(panel_data_pub,
     title    = 'Subcellular Proteome — Union Significant Proteins',
     footnote = ('AW-M-3 excluded  |  Union: corrected p (-log2) ≥ 3.3 in ≥1 condition  |  '
                 'Chromatin: Keep + Review only  |  Sorted by AW FC high → low'),
     out_pdf  = os.path.join(SCRIPT_DIR, 'Heatmap_Union_PubQuality_v2.pdf'),
     out_png  = os.path.join(SCRIPT_DIR, 'Heatmap_Union_PubQuality_v2.png'))

# Copy to subfolder
pub_dir = os.path.join(SCRIPT_DIR, 'Union_Heatmaps_PubQuality')
os.makedirs(pub_dir, exist_ok=True)
for ext in ['pdf', 'png']:
    shutil.copy2(os.path.join(SCRIPT_DIR, f'Heatmap_Union_PubQuality_v2.{ext}'),
                 os.path.join(pub_dir,     f'Heatmap_Union_PubQuality_v2.{ext}'))


# ══════════════════════════════════════════════════════════════════════════════
# HEATMAP 2 — Panel Union (computed fresh from raw reps)
# ══════════════════════════════════════════════════════════════════════════════
print('\n=== Heatmap 2: Panel Union v2 ===')
PREFIXES = {'Membrane':'Memb_', 'Cytosol':'Cyto_',
            'Chromatin':'Chrom_', 'Soluble Nuclear':'Nuc_'}

panel_data_panel = []
for disp, sheet, fmode in COMPARTMENTS_4:
    df = pd.read_excel(FILE, sheet_name=sheet)
    if fmode == 'keep_review' and 'Filter' in df.columns:
        df = df[df['Filter'].isin(['Keep','Review'])].reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    px   = PREFIXES[disp]
    n_c  = get_cols(df, px, NAIVE_SUFFIXES)
    i_c  = get_cols(df, px, INTOX_SUFFIXES)
    aw_c = get_cols(df, px, AW_SUFFIXES)
    pa_c = get_cols(df, px, PA_SUFFIXES)

    union_mask = pd.Series(False, index=df.index)
    for cc in [i_c, aw_c, pa_c]:
        if cc:
            fc   = mean_fc(df, n_c, cc)
            corr = corrected_p(df, n_c, cc)
            union_mask |= (corr > CORR_THRESH) & (fc.abs() > FC_THRESH)

    df_u = df[union_mask].reset_index(drop=True)
    n    = len(df_u)
    print(f'  {disp}: {n}')
    if n == 0: continue

    nc2  = get_cols(df_u, px, NAIVE_SUFFIXES)
    fi   = mean_fc(df_u, nc2, get_cols(df_u, px, INTOX_SUFFIXES)).values
    fa   = mean_fc(df_u, nc2, get_cols(df_u, px, AW_SUFFIXES)).values
    fp   = mean_fc(df_u, nc2, get_cols(df_u, px, PA_SUFFIXES)).values
    mat  = np.vstack([fi, fa, fp])
    order = np.argsort(mat[1])[::-1]
    panel_data_panel.append(dict(label=disp, n=n, mat=mat[:, order]))

draw(panel_data_panel,
     title    = 'Subcellular Proteome — Union Significant Proteins  (corrected p > 3.3, |FC| > 0.5)',
     footnote = ('AW-M-3 excluded  |  Union: corrected p (-log2) > 3.3 and |FC| > 0.5 in ≥1 condition  |  '
                 'Chromatin: Keep + Review only  |  Sorted by AW FC high → low'),
     out_pdf  = os.path.join(SCRIPT_DIR, 'Heatmap_Panel_Union_v2.pdf'),
     out_png  = os.path.join(SCRIPT_DIR, 'Heatmap_Panel_Union_v2.png'))


print('\nAll done.')
