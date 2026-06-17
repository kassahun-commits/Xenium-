"""
Paper figure — True translocation (Chromatin <-> Soluble Nuclear)
=================================================================
Two-panel publication figure:

  Panel A : TRUE-translocation UNION across all 3 conditions
            (proteins that truly translocate in >=1 condition).
            Same protein set shown across all conditions
            (Intox / AW / PA), each with Chromatin + SN columns.

  Panel B : AW-only TRUE translocation (NOT a union) — proteins that
            truly translocate during Acute Withdrawal.

To keep protein names readable, each panel shows only the strongest
translocators: the top-N "Into Chromatin" + top-N "Into SN"
(default N=25 per direction => 50 labelled proteins per panel).

True-translocation criteria (per condition):
  Into Chromatin : Sig(p_adj<P) AND interaction>0 AND SN_FC   < -FC
  Into SN        : Sig(p_adj<P) AND interaction<0 AND SN_FC > 0 AND Chrom_FC < -FC
  (interaction = mean(delta_chrom) - mean(delta_SN) per replicate; Welch t-test
   on the two delta vectors; BH FDR per condition.)

Standing rules (MEWS lab):
  * Workbook: "EDIT Excluding AWM3 ..." (AW-M-3 already excluded)
  * Chromatin filtered to Filter in {Keep, Review}; SN uses ALL proteins
  * Editable vector text (pdf.fonttype=42); source CSV exported alongside figure
  * No hardcoded paths — everything via CLI args

Usage:
  python3 make_paper_figure_trueTransloc_union_AW.py \
      --xlsx  "../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx" \
      --outdir current
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
CH_I_COLS  = ['Chrom_I-F-1', 'Chrom_I-F-2', 'Chrom_I-F-3', 'Chrom_I-M-1', 'Chrom_I-M-2']
SN_I_COLS  = ['Nuc_I-F-1',   'Nuc_I-F-2',   'Nuc_I-F-3',   'Nuc_I-M-1',   'Nuc_I-M-2']
CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1',  'Nuc_AW-F-2',  'Nuc_AW-M-1',  'Nuc_AW-M-2']
CH_PA_COLS = ['Chrom_PA-F-1', 'Chrom_PA-F-2', 'Chrom_PA-F-3', 'Chrom_PA-M-1', 'Chrom_PA-M-2']
SN_PA_COLS = ['Nuc_PA-F-1',  'Nuc_PA-F-2',  'Nuc_PA-F-3',  'Nuc_PA-M-1',  'Nuc_PA-M-2']

CONDITIONS  = [('Intox', CH_I_COLS, SN_I_COLS),
               ('AW',    CH_AW_COLS, SN_AW_COLS),
               ('PA',    CH_PA_COLS, SN_PA_COLS)]
COND_FULL = {'Intox': 'Intoxication', 'AW': 'Acute Withdrawal', 'PA': 'Protracted Abstinence'}


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


def build_stats(xlsx, p_thresh, fc_thresh):
    """Return df_all with per-condition Chrom_FC/SN_FC/Interaction/p_adj/Sig/Dir/True."""
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
        row = {'Gene': g}
        for nm, ch_cols, sn_cols in CONDITIONS:
            d_ch = get_vals(ch_df, g, ch_cols) - np.nanmean(ch_naive)
            d_sn = get_vals(sn_df, g, sn_cols) - np.nanmean(sn_naive)
            interaction = float(np.nanmean(d_ch) - np.nanmean(d_sn))
            a = d_ch[np.isfinite(d_ch)]
            b = d_sn[np.isfinite(d_sn)]
            p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
            row[f'Interaction_{nm}'] = round(interaction, 4) if np.isfinite(interaction) else np.nan
            row[f'p_{nm}'] = p
            row[f'Chrom_FC_{nm}'] = round(mean_fc(ch_df, g, ch_cols, CH_N_COLS), 4)
            row[f'SN_FC_{nm}'] = round(mean_fc(sn_df, g, sn_cols, SN_N_COLS), 4)
        rows.append(row)
    df = pd.DataFrame(rows)

    for nm, _, _ in CONDITIONS:
        valid = df[f'p_{nm}'].notna()
        df[f'p_adj_{nm}'] = np.nan
        df.loc[valid, f'p_adj_{nm}'] = bh_correction(df.loc[valid, f'p_{nm}'].values)
        df[f'Sig_{nm}'] = df[f'p_adj_{nm}'] < p_thresh
        sig = df[f'Sig_{nm}']
        inter = df[f'Interaction_{nm}']
        snfc = pd.to_numeric(df[f'SN_FC_{nm}'], errors='coerce')
        chfc = pd.to_numeric(df[f'Chrom_FC_{nm}'], errors='coerce')
        true_ch = sig & (inter > 0) & (snfc < -fc_thresh)
        true_sn = sig & (inter < 0) & (snfc > 0) & (chfc < -fc_thresh)
        df[f'TrueDir_{nm}'] = np.where(true_ch, 'Into_Chromatin',
                              np.where(true_sn, 'Into_SN', 'NS'))
    return df


def pick_top(df, sort_col, top_per_dir):
    """Top-N most positive (Into Chromatin) + top-N most negative (Into SN)."""
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


def dir_to_num(series):
    return series.map({'Into_Chromatin': 1, 'Into_SN': -1}).fillna(0).values


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--p-thresh', type=float, default=0.10)
    ap.add_argument('--fc-thresh', type=float, default=0.5)
    ap.add_argument('--top-per-dir', type=int, default=25,
                    help='proteins labelled per direction per panel (default 25 => 50/panel)')
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--date', default=_dt.date.today().isoformat())
    ap.add_argument('--file-tag', default='RatAlcoholProteome')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pthr_tag = str(args.p_thresh).replace('.', 'p')
    fc_tag = str(args.fc_thresh).replace('.', 'p')
    n_lab = 2 * args.top_per_dir
    base = (f'{args.file_tag}_TrueTranslocation_PaperFigure_UnionAndAW'
            f'_top{n_lab}_p{pthr_tag}_FC{fc_tag}_{args.date}')
    out_pdf = os.path.join(args.outdir, base + '.pdf')
    out_png = os.path.join(args.outdir, base + '.png')
    out_csv_a = os.path.join(args.outdir, base + '_PanelA_TrueUnion_source.csv')
    out_csv_b = os.path.join(args.outdir, base + '_PanelB_AW_source.csv')

    print('Building per-condition translocation stats...')
    df = build_stats(args.xlsx, args.p_thresh, args.fc_thresh)

    true_any = ((df['TrueDir_Intox'] != 'NS') |
                (df['TrueDir_AW'] != 'NS') |
                (df['TrueDir_PA'] != 'NS'))
    df_tu = df[true_any].copy()
    df_tu['sum_interaction'] = df_tu[['Interaction_Intox', 'Interaction_AW',
                                      'Interaction_PA']].fillna(0).sum(axis=1)
    df_aw = df[df['TrueDir_AW'] != 'NS'].copy()
    print(f'  True union (>=1 condition): {len(df_tu)}   AW true: {len(df_aw)}')

    panelA = pick_top(df_tu, 'sum_interaction', args.top_per_dir)
    panelB = pick_top(df_aw, 'Interaction_AW', args.top_per_dir)
    print(f'  Panel A labelled: {len(panelA)}   Panel B labelled: {len(panelB)}')

    # ── source CSVs ────────────────────────────────────────────────────────────
    a_cols = ['Gene', 'sum_interaction']
    for nm, _, _ in CONDITIONS:
        a_cols += [f'Chrom_FC_{nm}', f'SN_FC_{nm}', f'Interaction_{nm}',
                   f'p_adj_{nm}', f'TrueDir_{nm}']
    panelA[a_cols].to_csv(out_csv_a, index=False)
    b_cols = ['Gene', 'Chrom_FC_AW', 'SN_FC_AW', 'Interaction_AW',
              'p_adj_AW', 'TrueDir_AW']
    panelB[b_cols].to_csv(out_csv_b, index=False)
    print(f'  Saved: {out_csv_a}\n  Saved: {out_csv_b}')

    # ── layout (inches) ────────────────────────────────────────────────────────
    nA, nB = len(panelA), len(panelB)
    n_rows = max(nA, nB)
    ROW_H = 0.16
    HM_H = n_rows * ROW_H
    FC_COL_W = 0.26
    DIR_COL_W = 0.13
    DIR_GAP = 0.06
    LABEL_W = 1.45
    PANEL_GAP = 0.75
    LEFT_PAD = 0.15
    RIGHT_PAD = 0.35
    TOP_PAD = 1.05
    BOT_PAD = 1.45

    A_FC_W = 6 * FC_COL_W
    A_DIR_W = 3 * DIR_COL_W
    B_FC_W = 2 * FC_COL_W
    B_DIR_W = 1 * DIR_COL_W

    x_A_fc = LEFT_PAD + LABEL_W
    x_A_dir = x_A_fc + A_FC_W + DIR_GAP
    x_B_fc = x_A_dir + A_DIR_W + PANEL_GAP + LABEL_W
    x_B_dir = x_B_fc + B_FC_W + DIR_GAP

    FIG_W = x_B_dir + B_DIR_W + RIGHT_PAD
    FIG_H = BOT_PAD + HM_H + TOP_PAD
    y0 = BOT_PAD                       # heatmap bottom (inches)

    vmax = args.vmax
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

    def ax_in(x, y, w, h):
        return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])

    def draw_fc(ax_x, ax_w, sub, fc_cols, col_labels, group_spans):
        h = len(sub) * ROW_H
        yb = y0 + (HM_H - h)              # top-align shorter panel
        ax = ax_in(ax_x, yb, ax_w, h)
        mat = np.column_stack([pd.to_numeric(sub[c], errors='coerce').fillna(0).values
                               for c in fc_cols])
        im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                       interpolation='nearest')
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub['Gene'].tolist(), fontsize=5.0, fontfamily='Arial')
        ax.tick_params(axis='y', length=0, pad=2)
        ax.set_xticks(range(len(fc_cols)))
        ax.set_xticklabels(col_labels, fontsize=5.5, fontfamily='Arial')
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', length=0, pad=2)
        for c in range(1, len(fc_cols)):
            ax.axvline(c - 0.5, color='white', lw=0.6)
        for sp in ax.spines.values():
            sp.set_visible(False)
        # group headers above (axes-fraction y=1 is the TOP of the heatmap)
        for label, (c0, c1) in group_spans:
            xc = (c0 + c1) / 2.0
            ax.annotate(label, xy=(xc, 1), xycoords=('data', 'axes fraction'),
                        xytext=(0, 16), textcoords='offset points',
                        ha='center', va='bottom', fontsize=7,
                        fontweight='bold', fontfamily='Arial', color='#1A3A5C',
                        annotation_clip=False)
        return ax, im, yb, h

    def draw_dir(ax_x, ax_w, sub, dir_cols, col_labels, yb, h):
        ax = ax_in(ax_x, yb, ax_w, h)
        mat = np.column_stack([dir_to_num(sub[c]) for c in dir_cols])
        ax.imshow(mat, aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM,
                  interpolation='nearest')
        ax.set_yticks([])
        ax.set_xticks(range(len(dir_cols)))
        ax.set_xticklabels(col_labels, fontsize=5.0, fontfamily='Arial', rotation=90)
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', length=0, pad=2)
        for c in range(1, len(dir_cols)):
            ax.axvline(c - 0.5, color='white', lw=0.6)
        for sp in ax.spines.values():
            sp.set_visible(False)
        return ax

    # ── Panel A ────────────────────────────────────────────────────────────────
    a_fc_cols = ['Chrom_FC_Intox', 'SN_FC_Intox',
                 'Chrom_FC_AW', 'SN_FC_AW',
                 'Chrom_FC_PA', 'SN_FC_PA']
    a_col_labels = ['Chr', 'SN', 'Chr', 'SN', 'Chr', 'SN']
    a_groups = [('Intox', (0, 1)), ('AW', (2, 3)), ('PA', (4, 5))]
    axA, imA, ybA, hA = draw_fc(x_A_fc, A_FC_W, panelA, a_fc_cols,
                                a_col_labels, a_groups)
    draw_dir(x_A_dir, A_DIR_W, panelA,
             ['TrueDir_Intox', 'TrueDir_AW', 'TrueDir_PA'],
             ['I', 'AW', 'PA'], ybA, hA)
    fig.text((x_A_fc + A_FC_W / 2) / FIG_W, (ybA + hA + 0.46) / FIG_H,
             'A', ha='center', va='bottom', fontsize=15, fontweight='bold',
             fontfamily='Arial')
    fig.text((x_A_fc + A_FC_W / 2) / FIG_W, (ybA - 0.16) / FIG_H,
             f'True-translocation union\n(across conditions, n={len(df_tu)}; top {len(panelA)} shown)',
             ha='center', va='top', fontsize=6.5, fontfamily='Arial', color='#444444')

    # ── Panel B ────────────────────────────────────────────────────────────────
    b_fc_cols = ['Chrom_FC_AW', 'SN_FC_AW']
    axB, imB, ybB, hB = draw_fc(x_B_fc, B_FC_W, panelB, b_fc_cols,
                                ['Chr', 'SN'], [('AW', (0, 1))])
    draw_dir(x_B_dir, B_DIR_W, panelB, ['TrueDir_AW'], ['AW'], ybB, hB)
    fig.text((x_B_fc + B_FC_W / 2) / FIG_W, (ybB + hB + 0.46) / FIG_H,
             'B', ha='center', va='bottom', fontsize=15, fontweight='bold',
             fontfamily='Arial')
    fig.text((x_B_fc + B_FC_W / 2) / FIG_W, (ybB - 0.16) / FIG_H,
             f'AW true translocation\n(n={len(df_aw)}; top {len(panelB)} shown)',
             ha='center', va='top', fontsize=6.5, fontfamily='Arial', color='#444444')

    # ── colorbar ───────────────────────────────────────────────────────────────
    cbar_w = 2.6
    ax_cb = ax_in((FIG_W - cbar_w) / 2, 0.55, cbar_w, 0.13)
    cb = fig.colorbar(imA, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-vmax, -2, -1, 0, 1, 2, vmax])
    cb.ax.tick_params(labelsize=6)
    ax_cb.set_xlabel('Log2 fold change vs Naive  (within fraction)',
                     fontsize=7, fontfamily='Arial', color='#444444', labelpad=3)

    # ── direction legend ─────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(facecolor='#E8305A', edgecolor='none', label='Into Chromatin'),
        mpatches.Patch(facecolor='#2B7FD4', edgecolor='none', label='Into SN'),
        mpatches.Patch(facecolor='#E8E8E8', edgecolor='#AAAAAA', label='NS'),
    ]
    fig.legend(handles=legend_patches, loc='lower right',
               bbox_to_anchor=(0.99, 0.012), fontsize=6.5, framealpha=0.9,
               edgecolor='#CCCCCC', title='Direction strip', title_fontsize=6.5)

    # ── title + footnote ─────────────────────────────────────────────────────────
    fig.text(0.5, 1 - 0.30 / FIG_H,
             'True protein translocation: Chromatin ↔ Soluble Nuclear',
             ha='center', va='top', fontsize=10, fontweight='bold',
             fontfamily='Arial', color='#1A3A5C')
    fig.text(0.5, 0.012,
             f'True translocation: BH FDR p_adj < {args.p_thresh}, |FC| ≥ {args.fc_thresh}, '
             f'opposite-direction movement  |  AW-M-3 excluded  |  '
             f'Chromatin: Keep+Review  |  SN: all proteins  |  '
             f'rows sorted Into-Chromatin (top) → Into-SN (bottom)',
             ha='center', va='bottom', fontsize=4.8, fontfamily='Arial', color='#888888')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor='white')
    fig.savefig(out_png, dpi=300, facecolor='white')
    plt.close(fig)
    print(f'  Saved: {out_pdf}\n  Saved: {out_png}\nDone.')


if __name__ == '__main__':
    main()
