import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

FILE = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT  = 'AW_Combined_Panel.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_UP    = '#E8305A'
C_DOWN  = '#2B7FD4'
C_OTHER = '#C8C8C8'
C_NS    = '#EBEBEB'
C_GRAY  = '#AAAAAA'

SHEETS = [
    ('Membrane',        'Membrane',        'Memb_N-',  'Memb_I-',  'Memb_AW-',  'Memb_PA-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-',  'Cyto_I-',  'Cyto_AW-',  'Cyto_PA-',  'all'),
    ('Chromatin',       'Chromatin',       'Chrom_N-', 'Chrom_I-', 'Chrom_AW-', 'Chrom_PA-', 'keep_review'),
    ('Soluble nuclear', 'Soluble nuclear', 'Nuc_N-',   'Nuc_I-',   'Nuc_AW-',   'Nuc_PA-',   'keep_review'),
]

PATCH_IDS  = ['100', '010', '110', '001', '101', '011', '111']
AW_PATCHES = {'100', '110', '101', '111'}   # all regions that include circle A (AW)

# ── Stats helpers ──────────────────────────────────────────────────────────────
def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        p_vals.append(stats.ttest_ind(n, c, equal_var=False)[1]
                      if len(n) >= 2 and len(c) >= 2 else np.nan)
    p_vals    = pd.Series(p_vals, index=df.index)
    valid     = p_vals.notna()
    ranks     = p_vals[valid].rank(ascending=False)
    corrected = pd.Series(np.nan, index=df.index)
    corrected[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corrected

def get_sig_genes(df, naive_cols, cond_cols, direction):
    fc, corr = calc_stats(df, naive_cols, cond_cols)
    mask = fc.notna() & corr.notna() & (corr > CORR_THRESH)
    sig  = mask & (fc > FC_THRESH if direction == 'up' else fc < -FC_THRESH)
    return set(df.loc[sig, 'Gene symbol'].astype(str).tolist())

# ── Volcano panel ──────────────────────────────────────────────────────────────
def draw_volcano(ax, fc, corrected, highlight, title):
    mask = fc.notna() & corrected.notna()
    f, cr = fc[mask], corrected[mask]
    up   = (cr > CORR_THRESH) & (f >  FC_THRESH)
    down = (cr > CORR_THRESH) & (f < -FC_THRESH)
    ns   = ~up & ~down
    hi, lo, c_hi = (up, down, C_UP) if highlight == 'up' else (down, up, C_DOWN)

    ax.scatter(f[ns], cr[ns], s=8,  color=C_NS,    alpha=0.5,  linewidths=0, rasterized=True)
    ax.scatter(f[lo], cr[lo], s=14, color=C_OTHER, alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(f[hi], cr[hi], s=28, color=c_hi,    alpha=1.0,  linewidths=0, rasterized=True)

    ax.axhline(CORR_THRESH, color='#555555', linestyle='--', linewidth=1.4, alpha=0.6)
    ax.axvline( FC_THRESH,  color='#555555', linestyle='--', linewidth=1.4, alpha=0.6)
    ax.axvline(-FC_THRESH,  color='#555555', linestyle='--', linewidth=1.4, alpha=0.6)
    ax.axhline(0,           color='#999999', linestyle='-',  linewidth=0.6, alpha=0.4)
    ax.axvline(0,           color='#999999', linestyle='-',  linewidth=0.6, alpha=0.4)

    ax.set_xlabel('Fold Change',                fontsize=21, labelpad=8)
    ax.set_ylabel('Corrected p-value (−log₂)', fontsize=21, labelpad=8)
    ax.set_title(title, fontsize=20, fontweight='bold', pad=10)
    ax.tick_params(labelsize=16)
    ax.spines['top'].set_visible(False);  ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2); ax.spines['bottom'].set_linewidth(1.2)

    n_hi  = int(hi.sum())
    n_tot = int(mask.sum())
    direction = '↑' if highlight == 'up' else '↓'
    ax.text(0.98, 0.98, f'{direction}{n_hi}   n={n_tot}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=20, color='#333333', fontweight='bold')

    dir_label = 'Increased' if highlight == 'up' else 'Decreased'
    ax.legend(handles=[
        mpatches.Patch(color=c_hi,    label=f'{dir_label}  (n={n_hi})'),
        mpatches.Patch(color=C_OTHER, label='Other sig.'),
        mpatches.Patch(color=C_NS,    label='NS'),
    ], fontsize=15, loc='upper left', framealpha=0.85,
       edgecolor='#CCCCCC', frameon=True, borderpad=0.7)

# ── Venn panel ─────────────────────────────────────────────────────────────────
def subset_sizes(A, B, C):
    return (len(A-B-C), len(B-A-C), len((A&B)-C),
            len(C-A-B), len((A&C)-B), len((B&C)-A), len(A&B&C))

def make_layout_subsets(A, B, C):
    sA = max(1, len(A)) ** 0.5
    sB = max(1, len(B)) ** 0.5
    sC = max(1, len(C)) ** 0.5
    ov = min(sA, sB, sC) * 0.18
    return (max(0.05, sA-2*ov), max(0.05, sB-2*ov), ov,
            max(0.05, sC-2*ov), ov, ov, ov*0.4)

def draw_venn(ax, A, B, C, set_labels, highlight_color, label_size=20):
    actual = subset_sizes(A, B, C)
    layout = make_layout_subsets(A, B, C)
    v = venn3(subsets=layout, set_labels=set_labels, ax=ax)

    for patch_id in PATCH_IDS:
        patch = v.get_patch_by_id(patch_id)
        if patch:
            if patch_id in AW_PATCHES:
                patch.set_facecolor(highlight_color)
                patch.set_alpha(0.75)
            else:
                patch.set_facecolor(C_GRAY)
                patch.set_alpha(0.30)
            patch.set_edgecolor('none')

    for lbl, val in zip(v.subset_labels, actual):
        if lbl is not None:
            lbl.set_text(str(val))
            lbl.set_fontsize(label_size)
            lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    vc = venn3_circles(subsets=layout, ax=ax)
    for circle in (vc.patches if hasattr(vc, 'patches') else vc):
        circle.set_edgecolor('#333333')
        circle.set_linewidth(3.5)
        circle.set_facecolor('none')

    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(label_size + 2)
            lbl.set_fontfamily('Arial')
            lbl.set_fontweight('bold')
            lbl.set_color('#111111')

    # Zoom so circles fill well but labels have room to breathe
    xl, yl = ax.get_xlim(), ax.get_ylim()
    xc = (xl[0]+xl[1])/2;  yc = (yl[0]+yl[1])/2
    xr = (xl[1]-xl[0])*0.70; yr = (yl[1]-yl[0])*0.70
    ax.set_xlim(xc-xr, xc+xr);  ax.set_ylim(yc-yr, yc+yr)

# ── Build PDF ──────────────────────────────────────────────────────────────────
with PdfPages(OUT) as pdf:
    for label, excel_name, n_pat, i_pat, aw_pat, pa_pat, fmode in SHEETS:
        raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
        df  = raw[raw['Filter'].isin(['Keep','Review'])].reset_index(drop=True) \
              if fmode == 'keep_review' else raw.copy().reset_index(drop=True)

        naive_cols = [c for c in df.columns if n_pat  in str(c)]
        i_cols     = [c for c in df.columns if i_pat  in str(c)]
        aw_cols    = [c for c in df.columns if aw_pat in str(c)]
        pa_cols    = [c for c in df.columns if pa_pat in str(c)]

        fc, corrected = calc_stats(df, naive_cols, aw_cols)

        aw_up = get_sig_genes(df, naive_cols, aw_cols, 'up')
        i_up  = get_sig_genes(df, naive_cols, i_cols,  'up')
        pa_up = get_sig_genes(df, naive_cols, pa_cols, 'up')
        aw_dn = get_sig_genes(df, naive_cols, aw_cols, 'down')
        i_dn  = get_sig_genes(df, naive_cols, i_cols,  'down')
        pa_dn = get_sig_genes(df, naive_cols, pa_cols, 'down')

        # ── Figure layout ──────────────────────────────────────────────────────
        fig = plt.figure(figsize=(20, 17))
        fig.patch.set_facecolor('white')

        # Outer grid: 2 rows, VP row shorter than Venn row
        outer = gridspec.GridSpec(
            2, 1,
            figure=fig,
            height_ratios=[1.0, 1.4],
            hspace=0.22,
            top=0.93, bottom=0.03,
            left=0.05, right=0.98
        )

        # VP row: tighter horizontal spacing
        gs_vp = gridspec.GridSpecFromSubplotSpec(
            1, 2, subplot_spec=outer[0],
            wspace=0.22
        )

        # Venn row: wider horizontal spacing so labels don't collide
        gs_vn = gridspec.GridSpecFromSubplotSpec(
            1, 2, subplot_spec=outer[1],
            wspace=0.50
        )

        ax_vp_up   = fig.add_subplot(gs_vp[0, 0])
        ax_vp_dn   = fig.add_subplot(gs_vp[0, 1])
        ax_venn_up = fig.add_subplot(gs_vn[0, 0])
        ax_venn_dn = fig.add_subplot(gs_vn[0, 1])

        # ── Row 1: Volcano plots ───────────────────────────────────────────────
        draw_volcano(ax_vp_up, fc, corrected, 'up',
                     f'Acute Withdrawal vs Naïve — Increased')
        draw_volcano(ax_vp_dn, fc, corrected, 'down',
                     f'Acute Withdrawal vs Naïve — Decreased')

        # ── Row 2: Venn diagrams ───────────────────────────────────────────────
        draw_venn(ax_venn_up, aw_up, i_up, pa_up,
                  set_labels=(
                      f'AW Increased\n(n={len(aw_up)})',
                      f'Intox. Increased\n(n={len(i_up)})',
                      f'PA Increased\n(n={len(pa_up)})',
                  ),
                  highlight_color=C_UP)

        draw_venn(ax_venn_dn, aw_dn, i_dn, pa_dn,
                  set_labels=(
                      f'AW Decreased\n(n={len(aw_dn)})',
                      f'Intox. Decreased\n(n={len(i_dn)})',
                      f'PA Decreased\n(n={len(pa_dn)})',
                  ),
                  highlight_color=C_DOWN)

        # Remove axes borders from Venn panels
        for ax in [ax_venn_up, ax_venn_dn]:
            ax.set_axis_off()

        # ── Panel labels (A, B, C, D) ──────────────────────────────────────────
        for ax, letter in zip([ax_vp_up, ax_vp_dn, ax_venn_up, ax_venn_dn],
                               ['A', 'B', 'C', 'D']):
            ax.text(-0.06, 1.05, letter, transform=ax.transAxes,
                    fontsize=26, fontweight='bold', fontfamily='Arial',
                    va='top', ha='left', color='#111111')

        # ── Page title ────────────────────────────────────────────────────────
        fig.suptitle(
            f'{label}  —  Acute Withdrawal',
            fontsize=26, fontweight='bold', fontfamily='Arial',
            y=0.97, color='#111111'
        )

        pdf.savefig(fig, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f'{label}: done')

print(f'Saved: {OUT}')
