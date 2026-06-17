"""
AW true translocation — split two-block heatmap (paper figure)
==============================================================
Single condition (Acute Withdrawal), original orientation kept:
proteins are ROWS with gene names on the side. To make the figure WIDER
and shorter without flipping, the genes are split into two side-by-side
blocks:

  LEFT block  : Into Chromatin  (chromatin up, SN down)
  RIGHT block : Into Soluble Nuclear  (chromatin down, SN up)

Each block: gene names (left) + Chromatin and SN fold-change columns.
Default top-N per direction = 35  => 70 genes total, 35 rows per block.

True-translocation criteria (AW):
  Into Chromatin : Sig(p_adj<P) AND interaction>0 AND SN_FC   < -FC
  Into SN        : Sig(p_adj<P) AND interaction<0 AND SN_FC > 0 AND Chrom_FC < -FC

Standing rules (MEWS lab):
  * Workbook: "EDIT Excluding AWM3 ..." (AW-M-3 already excluded)
  * Chromatin filtered to Filter in {Keep, Review}; SN uses ALL proteins
  * Editable vector text (pdf.fonttype=42); source CSV exported alongside figure
  * No hardcoded paths — everything via CLI args

Usage:
  python3 make_AW_true_translocation_split_heatmap.py \
      --xlsx  "../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx" \
      --outdir current --top-per-dir 35
  # widen further with --col-w (inches per fold-change column)
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
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


# ── colormap ───────────────────────────────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'], N=512)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--p-thresh', type=float, default=0.10)
    ap.add_argument('--fc-thresh', type=float, default=0.5)
    ap.add_argument('--top-per-dir', type=int, default=35,
                    help='genes per direction / per block (default 35)')
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--col-w', type=float, default=0.60,
                    help='inches per fold-change column (raise to widen)')
    ap.add_argument('--row-h', type=float, default=0.17, help='inches per gene row')
    ap.add_argument('--date', default=_dt.date.today().isoformat())
    ap.add_argument('--file-tag', default='RatAlcoholProteome')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pthr_tag = str(args.p_thresh).replace('.', 'p')
    fc_tag = str(args.fc_thresh).replace('.', 'p')
    n_total = 2 * args.top_per_dir
    base = (f'{args.file_tag}_AWtrueTranslocation_SplitHeatmap'
            f'_top{n_total}_p{pthr_tag}_FC{fc_tag}_{args.date}')
    out_pdf = os.path.join(args.outdir, base + '.pdf')
    out_png = os.path.join(args.outdir, base + '.png')
    out_csv = os.path.join(args.outdir, base + '_source.csv')

    print('Building AW translocation stats...')
    df = build_aw_stats(args.xlsx, args.p_thresh, args.fc_thresh)
    df_aw = df[df['TrueDir_AW'] != 'NS'].copy()
    n_ch_all = int((df_aw['TrueDir_AW'] == 'Into_Chromatin').sum())
    n_sn_all = int((df_aw['TrueDir_AW'] == 'Into_SN').sum())
    print(f'  AW true translocation: {len(df_aw)}  (Into Chromatin={n_ch_all}, Into SN={n_sn_all})')

    into_ch = (df_aw[df_aw['TrueDir_AW'] == 'Into_Chromatin']
               .sort_values('Interaction_AW', ascending=False)
               .head(args.top_per_dir).reset_index(drop=True))
    into_sn = (df_aw[df_aw['TrueDir_AW'] == 'Into_SN']
               .sort_values('Interaction_AW', ascending=True)
               .head(args.top_per_dir).reset_index(drop=True))
    print(f'  Shown: Into Chromatin={len(into_ch)}, Into SN={len(into_sn)}')

    src = pd.concat([into_ch, into_sn])[['Gene', 'Chrom_FC_AW', 'SN_FC_AW',
                                         'Interaction_AW', 'p_adj_AW', 'TrueDir_AW']]
    src.to_csv(out_csv, index=False)
    print(f'  Saved: {out_csv}')

    # ── layout (inches) ────────────────────────────────────────────────────────
    ROW_H = args.row_h
    FC_COL_W = args.col_w
    LABEL_W = 1.30
    BLOCK_DATA_W = 2 * FC_COL_W
    BLOCK_GAP = 1.05
    LEFT_PAD = 0.20
    RIGHT_PAD = 0.35
    TOP_PAD = 1.05
    BOT_PAD = 1.35

    n_rows = max(len(into_ch), len(into_sn))
    HM_H = n_rows * ROW_H
    block_w = LABEL_W + BLOCK_DATA_W
    FIG_W = LEFT_PAD + block_w + BLOCK_GAP + block_w + RIGHT_PAD
    FIG_H = BOT_PAD + HM_H + TOP_PAD
    y0 = BOT_PAD

    vmax = args.vmax
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

    def ax_in(x, y, w, h):
        return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])

    def draw_block(x_label_left, sub, header, header_color):
        h = len(sub) * ROW_H
        yb = y0 + (HM_H - h)                       # top-align
        x_fc = x_label_left + LABEL_W
        ax = ax_in(x_fc, yb, BLOCK_DATA_W, h)
        mat = np.column_stack([
            pd.to_numeric(sub['Chrom_FC_AW'], errors='coerce').fillna(0).values,
            pd.to_numeric(sub['SN_FC_AW'], errors='coerce').fillna(0).values])
        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub['Gene'].tolist(), fontsize=6.0, fontfamily='Arial')
        ax.tick_params(axis='y', length=0, pad=2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Chromatin', 'Soluble\nNuclear'], fontsize=6.5,
                           fontfamily='Arial')
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', length=0, pad=3)
        ax.axvline(0.5, color='white', lw=1.2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # colored block header
        ax.annotate(header, xy=(0.5, 1), xycoords='axes fraction',
                    xytext=(0, 26), textcoords='offset points',
                    ha='center', va='bottom', fontsize=9, fontweight='bold',
                    fontfamily='Arial', color=header_color, annotation_clip=False)
        return im

    x_left = LEFT_PAD
    x_right = LEFT_PAD + block_w + BLOCK_GAP
    im = draw_block(x_left, into_ch, f'Into Chromatin  (n={len(into_ch)})', '#C01E42')
    draw_block(x_right, into_sn, f'Into Soluble Nuclear  (n={len(into_sn)})', '#1A5FA0')

    # Colorbar
    cbar_w = 2.6
    ax_cb = ax_in((FIG_W - cbar_w) / 2, 0.50, cbar_w, 0.13)
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-vmax, -2, -1, 0, 1, 2, vmax])
    cb.ax.tick_params(labelsize=6)
    ax_cb.set_xlabel('Log2 fold change vs Naive  (within fraction)',
                     fontsize=7, fontfamily='Arial', color='#444444', labelpad=3)

    # Title + footnote
    fig.text(0.5, 1 - 0.30 / FIG_H,
             f'Acute Withdrawal: true protein translocation (Chromatin vs Soluble Nuclear)  '
             f'— top {len(into_ch) + len(into_sn)} of {len(df_aw)}',
             ha='center', va='top', fontsize=10, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')
    fig.text(0.5, 0.012,
             f'True translocation: BH FDR p_adj < {args.p_thresh}, |FC| >= {args.fc_thresh}, '
             f'opposite-direction movement  |  AW-M-3 excluded  |  '
             f'Chromatin: Keep+Review  |  SN: all proteins  |  '
             f'each block sorted by translocation magnitude (strongest at top)',
             ha='center', va='bottom', fontsize=4.8, fontfamily='Arial', color='#888888')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor='white')
    fig.savefig(out_png, dpi=300, facecolor='white')
    plt.close(fig)
    print(f'  Saved: {out_pdf}\n  Saved: {out_png}\nDone.')


if __name__ == '__main__':
    main()
