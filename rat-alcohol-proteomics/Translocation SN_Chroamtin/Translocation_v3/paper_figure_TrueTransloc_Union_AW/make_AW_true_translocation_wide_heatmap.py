"""
AW true translocation — wide horizontal heatmap (paper figure)
==============================================================
Single panel: proteins that TRULY translocate between Chromatin and the
Soluble Nuclear fraction during Acute Withdrawal (AW).

Wide / short layout:
  * proteins are COLUMNS (gene names rotated at the bottom)
  * two heatmap rows: Chromatin (top) and Soluble Nuclear (bottom)
  * thin direction strip (red = Into Chromatin, blue = Into SN)
  * top-N most "Into Chromatin" + top-N most "Into SN"
    (default N=35 per direction => 70 genes)

True-translocation criteria (AW):
  Into Chromatin : Sig(p_adj<P) AND interaction>0 AND SN_FC   < -FC
  Into SN        : Sig(p_adj<P) AND interaction<0 AND SN_FC > 0 AND Chrom_FC < -FC
  (interaction = mean(delta_chrom) - mean(delta_SN) per replicate; Welch t-test;
   BH FDR across all common proteins.)

Standing rules (MEWS lab):
  * Workbook: "EDIT Excluding AWM3 ..." (AW-M-3 already excluded)
  * Chromatin filtered to Filter in {Keep, Review}; SN uses ALL proteins
  * Editable vector text (pdf.fonttype=42); source CSV exported alongside figure
  * No hardcoded paths — everything via CLI args

Usage:
  python3 make_AW_true_translocation_wide_heatmap.py \
      --xlsx  "../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx" \
      --outdir current --top-per-dir 35
"""

import os
import argparse
import datetime as _dt
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# ── Sample columns (AW-M-3 already absent from workbook) ───────────────────────
CH_N_COLS  = ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',   'Nuc_N-F-2',   'Nuc_N-F-3',   'Nuc_N-M-1',   'Nuc_N-M-2']
CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1',  'Nuc_AW-F-2',  'Nuc_AW-M-1',  'Nuc_AW-M-2']


def bh_correction(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    p_adj = pvals[order] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        p_adj[i] = min(p_adj[i], p_adj[i + 1])
    result = np.empty(n)
    result[order] = np.minimum(1.0, p_adj)
    return result


def get_vals(df, gene, cols):
    if gene in df.index:
        return pd.to_numeric(df.loc[gene, cols], errors='coerce').values
    return np.full(len(cols), np.nan)


def mean_fc(df, gene, cond_cols, naive_cols):
    return float(np.nanmean(get_vals(df, gene, cond_cols)) -
                 np.nanmean(get_vals(df, gene, naive_cols)))


def build_aw_stats(xlsx, p_thresh, fc_thresh):
    ch_raw = pd.read_excel(xlsx, sheet_name='Chromatin')
    ch_df = ch_raw[ch_raw['Filter'].isin(['Keep', 'Review'])].copy()
    sn_df = pd.read_excel(xlsx, sheet_name='Soluble nuclear').copy()
    for d in (ch_df, sn_df):
        d['Gene symbol'] = d['Gene symbol'].astype(str).str.strip()
    ch_df = ch_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
    sn_df = sn_df.drop_duplicates('Gene symbol').set_index('Gene symbol')
    common = sorted(set(ch_df.index) & set(sn_df.index))

    rows = []
    for g in common:
        ch_naive = get_vals(ch_df, g, CH_N_COLS)
        sn_naive = get_vals(sn_df, g, SN_N_COLS)
        d_ch = get_vals(ch_df, g, CH_AW_COLS) - np.nanmean(ch_naive)
        d_sn = get_vals(sn_df, g, SN_AW_COLS) - np.nanmean(sn_naive)
        interaction = float(np.nanmean(d_ch) - np.nanmean(d_sn))
        a = d_ch[np.isfinite(d_ch)]
        b = d_sn[np.isfinite(d_sn)]
        p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
        rows.append({'Gene': g,
                     'Interaction_AW': round(interaction, 4) if np.isfinite(interaction) else np.nan,
                     'p_AW': p,
                     'Chrom_FC_AW': round(mean_fc(ch_df, g, CH_AW_COLS, CH_N_COLS), 4),
                     'SN_FC_AW': round(mean_fc(sn_df, g, SN_AW_COLS, SN_N_COLS), 4)})
    df = pd.DataFrame(rows)

    valid = df['p_AW'].notna()
    df['p_adj_AW'] = np.nan
    df.loc[valid, 'p_adj_AW'] = bh_correction(df.loc[valid, 'p_AW'].values)
    sig = df['p_adj_AW'] < p_thresh
    inter = df['Interaction_AW']
    snfc = pd.to_numeric(df['SN_FC_AW'], errors='coerce')
    chfc = pd.to_numeric(df['Chrom_FC_AW'], errors='coerce')
    true_ch = sig & (inter > 0) & (snfc < -fc_thresh)
    true_sn = sig & (inter < 0) & (snfc > 0) & (chfc < -fc_thresh)
    df['TrueDir_AW'] = np.where(true_ch, 'Into_Chromatin',
                       np.where(true_sn, 'Into_SN', 'NS'))
    return df


def pick_top(df, sort_col, top_per_dir):
    s = df.copy()
    s[sort_col] = pd.to_numeric(s[sort_col], errors='coerce').fillna(0)
    into_ch = s.sort_values(sort_col, ascending=False).head(top_per_dir)
    into_sn = s.sort_values(sort_col, ascending=True).head(top_per_dir)
    out = pd.concat([into_ch, into_sn]).drop_duplicates('Gene')
    return out.sort_values(sort_col, ascending=False).reset_index(drop=True)


# ── colormaps ──────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'], N=512)
ANNOT_CMAP = ListedColormap(['#2B7FD4', '#E8E8E8', '#E8305A'])   # SN / NS / Chromatin
ANNOT_NORM = BoundaryNorm([-1.5, -0.5, 0.5, 1.5], ANNOT_CMAP.N)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--p-thresh', type=float, default=0.10)
    ap.add_argument('--fc-thresh', type=float, default=0.5)
    ap.add_argument('--top-per-dir', type=int, default=35,
                    help='genes per direction (default 35 => 70 total)')
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--col-w', type=float, default=0.17, help='inches per gene column')
    ap.add_argument('--date', default=_dt.date.today().isoformat())
    ap.add_argument('--file-tag', default='RatAlcoholProteome')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pthr_tag = str(args.p_thresh).replace('.', 'p')
    fc_tag = str(args.fc_thresh).replace('.', 'p')
    n_lab = 2 * args.top_per_dir
    base = (f'{args.file_tag}_AWtrueTranslocation_WideHeatmap'
            f'_top{n_lab}_p{pthr_tag}_FC{fc_tag}_{args.date}')
    out_pdf = os.path.join(args.outdir, base + '.pdf')
    out_png = os.path.join(args.outdir, base + '.png')
    out_csv = os.path.join(args.outdir, base + '_source.csv')

    print('Building AW translocation stats...')
    df = build_aw_stats(args.xlsx, args.p_thresh, args.fc_thresh)
    df_aw = df[df['TrueDir_AW'] != 'NS'].copy()
    n_ch = int((df_aw['TrueDir_AW'] == 'Into_Chromatin').sum())
    n_sn = int((df_aw['TrueDir_AW'] == 'Into_SN').sum())
    print(f'  AW true translocation: {len(df_aw)}  (Into Chromatin={n_ch}, Into SN={n_sn})')

    sub = pick_top(df_aw, 'Interaction_AW', args.top_per_dir)
    n = len(sub)
    print(f'  Shown: {n}')
    sub[['Gene', 'Chrom_FC_AW', 'SN_FC_AW', 'Interaction_AW',
         'p_adj_AW', 'TrueDir_AW']].to_csv(out_csv, index=False)
    print(f'  Saved: {out_csv}')

    # ── layout (inches) ────────────────────────────────────────────────────────
    COL_W = args.col_w
    HM_W = n * COL_W
    ROW_FC = 0.46
    ROW_DIR = 0.17
    ROW_GAP = 0.05
    LEFT = 1.70
    RIGHT = 1.95
    TOP = 0.85
    LABEL_BAND = 1.20      # rotated gene names below the heatmap
    CBAR_BAND = 0.95

    FIG_W = LEFT + HM_W + RIGHT
    FIG_H = TOP + ROW_FC + ROW_FC + ROW_GAP + ROW_DIR + LABEL_BAND + CBAR_BAND

    vmax = args.vmax
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

    def ax_in(x, y, w, h):
        return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])

    x_hm = LEFT
    y_top = FIG_H - TOP
    y_ch = y_top - ROW_FC                      # chromatin row bottom
    y_sn = y_ch - ROW_FC                        # SN row bottom
    y_dir = y_sn - ROW_GAP - ROW_DIR           # direction strip bottom

    chrom_mat = pd.to_numeric(sub['Chrom_FC_AW'], errors='coerce').fillna(0).values.reshape(1, -1)
    sn_mat = pd.to_numeric(sub['SN_FC_AW'], errors='coerce').fillna(0).values.reshape(1, -1)
    dir_mat = sub['TrueDir_AW'].map({'Into_Chromatin': 1, 'Into_SN': -1}).fillna(0).values.reshape(1, -1)

    # Chromatin row
    ax_ch = ax_in(x_hm, y_ch, HM_W, ROW_FC)
    im = ax_ch.imshow(chrom_mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_ch.set_xticks([]); ax_ch.set_yticks([])
    for sp in ax_ch.spines.values(): sp.set_visible(False)

    # SN row
    ax_sn = ax_in(x_hm, y_sn, HM_W, ROW_FC)
    ax_sn.imshow(sn_mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax_sn.set_xticks([]); ax_sn.set_yticks([])
    for sp in ax_sn.spines.values(): sp.set_visible(False)

    # white line between the two fraction rows
    fig.add_artist(plt.Line2D([x_hm / FIG_W, (x_hm + HM_W) / FIG_W],
                              [y_ch / FIG_H, y_ch / FIG_H],
                              transform=fig.transFigure, color='white', lw=1.2, zorder=5))

    # Direction strip
    ax_dir = ax_in(x_hm, y_dir, HM_W, ROW_DIR)
    ax_dir.imshow(dir_mat, aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM, interpolation='nearest')
    ax_dir.set_yticks([])
    ax_dir.set_xticks(range(n))
    ax_dir.set_xticklabels(sub['Gene'].tolist(), fontsize=5.2, fontfamily='Arial',
                           rotation=90, va='top')
    ax_dir.tick_params(axis='x', length=0, pad=2)
    for sp in ax_dir.spines.values(): sp.set_visible(False)

    # boundary between the Into-Chromatin block and the Into-SN block
    n_into_ch = int((sub['TrueDir_AW'] == 'Into_Chromatin').sum())
    if 0 < n_into_ch < n:
        xb = (x_hm + n_into_ch * COL_W) / FIG_W
        fig.add_artist(plt.Line2D([xb, xb],
                                  [y_dir / FIG_H, y_top / FIG_H],
                                  transform=fig.transFigure, color='#333333',
                                  lw=1.0, ls=(0, (4, 2)), zorder=6))

    # Row labels (left)
    for yb, lab in [(y_ch, 'Chromatin'), (y_sn, 'Soluble\nNuclear')]:
        fig.text((x_hm - 0.12) / FIG_W, (yb + ROW_FC / 2) / FIG_H, lab,
                 ha='right', va='center', fontsize=9, fontweight='bold',
                 fontfamily='Arial', color='#111111')
    fig.text((x_hm - 0.12) / FIG_W, (y_dir + ROW_DIR / 2) / FIG_H, 'Direction',
             ha='right', va='center', fontsize=6.5, fontfamily='Arial',
             color='#555555', style='italic')

    # Block headers above the heatmap
    fig.text((x_hm + n_into_ch * COL_W / 2) / FIG_W, (y_top + 0.06) / FIG_H,
             f'Into Chromatin  (n={n_into_ch})', ha='center', va='bottom',
             fontsize=8, fontweight='bold', fontfamily='Arial', color='#C01E42')
    fig.text((x_hm + (n_into_ch + (n - n_into_ch) / 2) * COL_W) / FIG_W, (y_top + 0.06) / FIG_H,
             f'Into Soluble Nuclear  (n={n - n_into_ch})', ha='center', va='bottom',
             fontsize=8, fontweight='bold', fontfamily='Arial', color='#1A5FA0')

    # Colorbar
    cbar_w = 2.6
    ax_cb = ax_in((FIG_W - cbar_w) / 2, 0.42, cbar_w, 0.13)
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-vmax, -2, -1, 0, 1, 2, vmax])
    cb.ax.tick_params(labelsize=6)
    ax_cb.set_xlabel('Log2 fold change vs Naive  (within fraction)',
                     fontsize=7, fontfamily='Arial', color='#444444', labelpad=3)

    # Direction legend
    legend_patches = [
        mpatches.Patch(facecolor='#E8305A', edgecolor='none', label='Into Chromatin'),
        mpatches.Patch(facecolor='#2B7FD4', edgecolor='none', label='Into SN'),
    ]
    fig.legend(handles=legend_patches, loc='lower right',
               bbox_to_anchor=(0.995, 0.02), fontsize=6.5, framealpha=0.9,
               edgecolor='#CCCCCC', title='Direction strip', title_fontsize=6.5)

    # Title + footnote
    fig.text(0.5, 1 - 0.30 / FIG_H,
             f'Acute Withdrawal: true protein translocation (Chromatin vs Soluble Nuclear)  '
             f'— top {n} of {len(df_aw)}',
             ha='center', va='top', fontsize=10, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')
    fig.text(0.5, 0.012,
             f'True translocation: BH FDR p_adj < {args.p_thresh}, |FC| >= {args.fc_thresh}, '
             f'opposite-direction movement  |  AW-M-3 excluded  |  '
             f'Chromatin: Keep+Review  |  SN: all proteins  |  '
             f'columns sorted Into-Chromatin (left) -> Into-SN (right)',
             ha='center', va='bottom', fontsize=4.8, fontfamily='Arial', color='#888888')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor='white')
    fig.savefig(out_png, dpi=300, facecolor='white')
    plt.close(fig)
    print(f'  Saved: {out_pdf}\n  Saved: {out_png}\nDone.')


if __name__ == '__main__':
    main()
