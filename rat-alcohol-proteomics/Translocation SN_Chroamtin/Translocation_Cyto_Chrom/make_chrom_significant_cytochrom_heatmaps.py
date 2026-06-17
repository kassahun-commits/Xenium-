#!/usr/bin/env python3
"""
Significant proteins (cytosol OR chromatin) shown across both compartments.
=============================================================================
This is the *exploration* companion to the strict translocation figure.  Here we
do NOT impose any cytosol-direction requirement and NO delta-delta gate.  A
protein is selected for a condition when it is significantly DIFFERENT (up OR
down) vs Naive in EITHER compartment -- the per-group UNION of the cytosol-
significant and chromatin-significant sets -- and we then display both its
Chromatin and Cytosol fold changes side by side.  From this you can eyeball which
proteins went up on chromatin while falling in cytosol (the true-translocation
candidates) and build the final translocation figure yourself.

SELECTION (criterion A in each compartment; default --select-mode union):
    chrom_sig = Chrom_Corrected >= CORR_THRESH  AND  |Chrom_FC| >= FC_THRESH
    cyto_sig  = Cyto_Corrected  >= CORR_THRESH  AND  |Cyto_FC|  >= FC_THRESH
    selected  = chrom_sig OR cyto_sig        (union; direction NOT filtered)
  --select-mode also supports: both (AND/intersection), chrom, cyto.

Two products are written:
  1. UNION-OVER-CONDITIONS heatmap -- every protein selected in >=1 condition,
                       displayed across all three condition bands (Chromatin row
                       over Cytosol row in each band).
  2. PER-CONDITION heatmaps (NOT union over conditions) -- one figure per
                       condition showing every protein selected in THAT condition,
                       with its Chromatin and Cytosol fold change.

Each protein's significance source per condition (chrom-only / cyto-only / both)
is shown as a categorical annotation strip.

Reported (CSV only, never used as a gate): Cyto_FC, Cyto_Corrected, the
delta-delta interaction score + Welch p / BH p_adj, and a convenience
'transloc_candidate' flag (chrom up & cyto down, or chrom down & cyto up).

Fold change = log2 vs Naive from the sheet's 'Fold change' columns (verified ==
mean(cond)-mean(naive) of the AW-M-3-excluded per-sample columns).  'Corrected'
is the sheet significance score (~4.0 at p~=0.05; 3.3 ~ p<=0.07).

Filtering: CHROMATIN restricted to Filter in {Keep, Review}.  Cytosol is NOT
filtered.  AW-M-3 already excluded from the workbook's fraction sheets.

MEWS rules: no hardcoded paths (CLI), editable vector text (pdf.fonttype=42),
source-data CSV alongside figures.  Outputs per product: PNG (300 dpi), editable
PDF, and CSVs (full per-protein stats + the exact plotted matrix).
"""
from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ---- per-sample log2 abundance columns (AW-M-3 already excluded) ------------
CH_SAMPLES = {
    'N':  ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2'],
    'I':  ['Chrom_I-F-1', 'Chrom_I-F-2', 'Chrom_I-F-3', 'Chrom_I-M-1', 'Chrom_I-M-2'],
    'AW': ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2'],
    'PA': ['Chrom_PA-F-1', 'Chrom_PA-F-2', 'Chrom_PA-F-3', 'Chrom_PA-M-1', 'Chrom_PA-M-2'],
}
CY_SAMPLES = {
    'N':  ['Cyto_N-F-1', 'Cyto_N-F-2', 'Cyto_N-F-3', 'Cyto_N-M-1', 'Cyto_N-M-2'],
    'I':  ['Cyto_I-F-1', 'Cyto_I-F-2', 'Cyto_I-F-3', 'Cyto_I-M-1', 'Cyto_I-M-2'],
    'AW': ['Cyto_AW-F-1', 'Cyto_AW-F-2', 'Cyto_AW-M-1', 'Cyto_AW-M-2'],
    'PA': ['Cyto_PA-F-1', 'Cyto_PA-F-2', 'Cyto_PA-F-3', 'Cyto_PA-M-1', 'Cyto_PA-M-2'],
}
SHEET_FC   = {'I': 'Fold change',   'AW': 'Fold change.1', 'PA': 'Fold change.2'}
SHEET_CORR = {'I': 'Corrected',     'AW': 'Corrected.1',   'PA': 'Corrected.2'}

CONDITIONS  = ['I', 'AW', 'PA']
COND_LABELS = {'I': 'Intoxication', 'AW': 'Acute Withdrawal', 'PA': 'Protracted Abstinence'}
COND_SHORT  = {'I': 'Intoxication', 'AW': 'AcuteWithdrawal', 'PA': 'ProtractedAbstinence'}

CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'], N=512)


def bh_correction(pvals):
    """Benjamini-Hochberg FDR; NaNs ignored (returned as NaN)."""
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    adj = pv[order] * n / (np.arange(n) + 1)
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    res = np.empty(n)
    res[order] = np.minimum(1.0, adj)
    out[ok] = res
    return out


def vals(df, cols):
    """Numeric matrix (rows = proteins) for the given sample columns."""
    return df[cols].apply(pd.to_numeric, errors='coerce').values


def build_stats(ch, cy, common):
    """Per-condition long table: Chrom/Cyto FC + Corrected, ddelta, p, p_adj, sig flag."""
    ch_naive_mean = np.nanmean(vals(ch, CH_SAMPLES['N']), axis=1)
    cy_naive_mean = np.nanmean(vals(cy, CY_SAMPLES['N']), axis=1)

    rows = []
    for cond in CONDITIONS:
        d_ch = vals(ch, CH_SAMPLES[cond]) - ch_naive_mean[:, None]
        d_cy = vals(cy, CY_SAMPLES[cond]) - cy_naive_mean[:, None]
        interaction = np.nanmean(d_ch, axis=1) - np.nanmean(d_cy, axis=1)

        pvec = np.full(len(common), np.nan)
        for i in range(len(common)):
            a = d_ch[i][np.isfinite(d_ch[i])]
            b = d_cy[i][np.isfinite(d_cy[i])]
            if a.size >= 2 and b.size >= 2:
                pvec[i] = stats.ttest_ind(a, b, equal_var=False).pvalue
        padj = bh_correction(pvec)

        chrom_fc   = pd.to_numeric(ch[SHEET_FC[cond]], errors='coerce').values
        cyto_fc    = pd.to_numeric(cy[SHEET_FC[cond]], errors='coerce').values
        chrom_corr = pd.to_numeric(ch[SHEET_CORR[cond]], errors='coerce').values
        cyto_corr  = pd.to_numeric(cy[SHEET_CORR[cond]], errors='coerce').values

        for i, acc in enumerate(common):
            rows.append({
                'Accession': acc,
                'Gene symbol': ch['Gene symbol'].iloc[i],
                'condition': COND_LABELS[cond],
                'cond_code': cond,
                'Chrom_FC': chrom_fc[i],
                'Chrom_Corrected': chrom_corr[i],
                'Cyto_FC': cyto_fc[i],
                'Cyto_Corrected': cyto_corr[i],
                'interaction_ddelta': interaction[i],
                'welch_p': pvec[i],
                'welch_p_adj': padj[i],
            })
    df = pd.DataFrame(rows)
    return df


def draw_band(fig, FIG_W, FIG_H, x0_in, y_bot_in, hm_w, frac_row_h,
              chrom_row, cyto_row, norm, label_right=True):
    """Draw one Chromatin-over-Cytosol band; returns the chromatin imshow handle."""
    COMP_BG = {'Chromatin': '#FBF9F0', 'Cytosol': '#F0F6FB'}
    ax_ch = fig.add_axes([x0_in / FIG_W, (y_bot_in + frac_row_h) / FIG_H,
                          hm_w / FIG_W, frac_row_h / FIG_H], facecolor=COMP_BG['Chromatin'])
    im = ax_ch.imshow(chrom_row.reshape(1, -1), aspect='auto', cmap=CMAP, norm=norm,
                      interpolation='nearest')
    ax_ch.set_xticks([]); ax_ch.set_yticks([])
    for sp in ax_ch.spines.values():
        sp.set_visible(False)

    ax_cy = fig.add_axes([x0_in / FIG_W, y_bot_in / FIG_H,
                          hm_w / FIG_W, frac_row_h / FIG_H], facecolor=COMP_BG['Cytosol'])
    ax_cy.imshow(cyto_row.reshape(1, -1), aspect='auto', cmap=CMAP, norm=norm,
                 interpolation='nearest')
    ax_cy.set_xticks([]); ax_cy.set_yticks([])
    for sp in ax_cy.spines.values():
        sp.set_visible(False)

    if label_right:
        fig.text((x0_in + hm_w + 0.10) / FIG_W,
                 (y_bot_in + frac_row_h + frac_row_h / 2) / FIG_H, 'Chromatin',
                 ha='left', va='center', fontsize=7.5, color='#111111')
        fig.text((x0_in + hm_w + 0.10) / FIG_W, (y_bot_in + frac_row_h / 2) / FIG_H,
                 'Cytosol', ha='left', va='center', fontsize=7.5, color='#111111')
    sep_y = (y_bot_in + frac_row_h) / FIG_H
    fig.add_artist(plt.Line2D([x0_in / FIG_W, (x0_in + hm_w) / FIG_W], [sep_y, sep_y],
                              transform=fig.transFigure, color='#CCCCCC', linewidth=0.8, zorder=5))
    return im


def hm_width(n_prot):
    """Adaptive heatmap width in inches (dense, no per-protein labels)."""
    return float(np.clip(n_prot * 0.045, 4.0, 18.0))


def write_plotmatrix(path, order, gene_of, sub, conds):
    out_rows = []
    look = {(r['Accession'], r['cond_code']): r for _, r in sub.iterrows()}
    for acc in order:
        for cond in conds:
            r = look.get((acc, cond))
            out_rows.append({
                'Accession': acc, 'Gene symbol': gene_of.get(acc, ''),
                'condition': COND_LABELS[cond],
                'Chrom_FC': r['Chrom_FC'] if r is not None else np.nan,
                'Chrom_Corrected': r['Chrom_Corrected'] if r is not None else np.nan,
                'Cyto_FC': r['Cyto_FC'] if r is not None else np.nan,
                'Cyto_Corrected': r['Cyto_Corrected'] if r is not None else np.nan,
                'chrom_sig': bool(r['chrom_sig']) if r is not None else False,
                'cyto_sig': bool(r['cyto_sig']) if r is not None else False,
                'selected': bool(r['selected']) if r is not None else False,
                'sig_source': r['sig_source_label'] if r is not None else 'none',
                'transloc_candidate': bool(r['transloc_candidate']) if r is not None else False,
            })
    pd.DataFrame(out_rows).to_csv(path, index=False)


SRC_COLORS = ['#E8E8E8', '#C0504D', '#4F81BD', '#7A4FB5']   # none, chrom, cyto, both
SRC_NAMES  = {1: 'chromatin only', 2: 'cytosol only', 3: 'both compartments'}


def plot_union(stats_df, args, outdir, base_tag):
    """Union over conditions: proteins selected in >=1 condition, three bands."""
    union = sorted(set(stats_df.loc[stats_df['selected'], 'Accession']))
    if not union:
        logging.warning('Union: no selected proteins; skipped.')
        return
    sub = stats_df[stats_df['Accession'].isin(union)].copy()
    gene_of = dict(zip(sub['Accession'], sub['Gene symbol']))

    def mat(field):
        m = sub.pivot(index='Accession', columns='cond_code', values=field).reindex(union)
        return m[CONDITIONS]

    chrom_mat = mat('Chrom_FC')
    cyto_mat  = mat('Cyto_FC')
    # significance-source matrix: 0 not-selected, 1 chrom-only, 2 cyto-only, 3 both
    src_mat = sub.copy()
    src_mat['src_sel'] = np.where(src_mat['selected'], src_mat['sig_source'], 0)
    src_mat = src_mat.pivot(index='Accession', columns='cond_code',
                            values='src_sel').reindex(union)[CONDITIONS].fillna(0).astype(int)

    order = chrom_mat.fillna(0).mean(axis=1).sort_values(ascending=False).index.tolist()
    chrom_mat = chrom_mat.reindex(order); cyto_mat = cyto_mat.reindex(order)
    src_mat = src_mat.reindex(order)
    n_prot = len(order)

    base = f'{base_tag}_Union'
    write_plotmatrix(outdir / f'{base}_plotmatrix.csv', order, gene_of, sub, CONDITIONS)

    norm = TwoSlopeNorm(vmin=-args.vmax, vcenter=0, vmax=args.vmax)
    ANNOT_CMAP = ListedColormap(SRC_COLORS)
    ANNOT_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], ANNOT_CMAP.N)

    FRAC_ROW_H = 0.50; BAND_H = 2 * FRAC_ROW_H
    ANNOT_ROW_H = 0.16; ANNOT_H = 3 * ANNOT_ROW_H
    ANNOT_GAP = 0.20; SECTION_GAP = 0.30; N_COND = 3
    LEFT_IN = 2.70; RIGHT_IN = 2.90; TOP_PAD = 0.70; BOT_PAD = 1.25
    CBAR_W_IN = 2.8; CBAR_H_IN = 0.14; CBAR_BOT = 0.46

    hm_w = hm_width(n_prot)
    FIG_W = LEFT_IN + hm_w + RIGHT_IN
    FIG_H = (TOP_PAD + ANNOT_H + ANNOT_GAP
             + N_COND * BAND_H + (N_COND - 1) * SECTION_GAP + BOT_PAD)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    y_cursor = FIG_H - TOP_PAD

    annot = np.vstack([src_mat[c].astype(int).values for c in CONDITIONS])
    y_annot_bot = y_cursor - ANNOT_H
    ax_an = fig.add_axes([LEFT_IN / FIG_W, y_annot_bot / FIG_H, hm_w / FIG_W, ANNOT_H / FIG_H])
    ax_an.imshow(annot, aspect='auto', cmap=ANNOT_CMAP, norm=ANNOT_NORM, interpolation='nearest')
    for r in range(1, N_COND):
        ax_an.axhline(r - 0.5, color='white', lw=0.8, zorder=3)
    ax_an.set_xticks([]); ax_an.set_yticks([])
    for sp in ax_an.spines.values():
        sp.set_visible(False)
    for r, cond in enumerate(CONDITIONS):
        yf = (y_annot_bot + ANNOT_H - (r + 0.5) * ANNOT_ROW_H) / FIG_H
        fig.text((LEFT_IN + hm_w + 0.10) / FIG_W, yf, COND_LABELS[cond],
                 ha='left', va='center', fontsize=6, color='#333333')
    fig.text((LEFT_IN - 0.12) / FIG_W, (y_annot_bot + ANNOT_H / 2) / FIG_H,
             'significance\nsource', ha='right', va='center',
             fontsize=7, color='#555555', style='italic')
    y_cursor = y_annot_bot - ANNOT_GAP

    im = None
    for cond in CONDITIONS:
        y_band_bot = y_cursor - BAND_H
        im = draw_band(fig, FIG_W, FIG_H, LEFT_IN, y_band_bot, hm_w, FRAC_ROW_H,
                       chrom_mat[cond].fillna(0).values, cyto_mat[cond].fillna(0).values, norm)
        fig.text((LEFT_IN - 0.12) / FIG_W, (y_band_bot + BAND_H / 2) / FIG_H,
                 COND_LABELS[cond].replace(' ', '\n'), ha='right', va='center',
                 fontsize=10, fontweight='bold', color='#111111')
        y_cursor = y_band_bot - SECTION_GAP

    if n_prot <= 80:
        bottom_y = (y_cursor + SECTION_GAP - 0.02) / FIG_H
        for pi, acc in enumerate(order):
            xf = (LEFT_IN + (pi + 0.5) / n_prot * hm_w) / FIG_W
            fig.text(xf, bottom_y, gene_of.get(acc, acc), ha='right', va='top',
                     rotation=90, fontsize=5.5, color='#222222', style='italic')

    legend_patches = [mpatches.Patch(color=SRC_COLORS[1], label='chromatin only'),
                      mpatches.Patch(color=SRC_COLORS[2], label='cytosol only'),
                      mpatches.Patch(color=SRC_COLORS[3], label='both compartments'),
                      mpatches.Patch(color=SRC_COLORS[0], label='not selected',
                                     linewidth=0.5, edgecolor='#AAAAAA')]
    fig.legend(handles=legend_patches, loc='lower right',
               bbox_to_anchor=(0.99, CBAR_BOT / FIG_H + 0.05),
               fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
               title='Significant in…', title_fontsize=7)

    cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
    ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H, CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3]); cb.ax.tick_params(labelsize=7)
    ax_cb.set_xlabel('log2 fold change vs Naive', fontsize=8, color='#444444', labelpad=4)

    fig.text(0.5, (FIG_H - TOP_PAD * 0.35) / FIG_H,
             f'Significant proteins ({args.select_mode}: cytosol/chromatin) by condition  '
             f'(union over conditions n = {n_prot})',
             ha='center', va='top', fontsize=11, fontweight='bold', color='#1A3A5C')
    fig.text(0.5, 0.008,
             f'Selection per condition: |FC|≥{args.fc_thresh:g} & Corrected≥{args.corr_thresh:g} '
             f'in {("chromatin OR cytosol" if args.select_mode=="union" else args.select_mode)} '
             f'(up OR down).  No ΔΔ filter.  '
             f'Chromatin filtered to {"/".join(args.keep_filter)}; cytosol unfiltered; AW-M-3 excluded.  '
             f'Proteins sorted by mean Chrom_FC; shown if selected in ≥1 condition.',
             ha='center', va='bottom', fontsize=5.5, color='#888888')

    png = outdir / f'{base}.png'; pdf = outdir / f'{base}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logging.info('UNION: wrote %s | %s (n=%d)', png.name, pdf.name, n_prot)


def plot_per_condition(stats_df, cond, args, outdir, base_tag):
    """One condition: every selected protein (sig in chrom OR cyto), both rows."""
    sub = stats_df[(stats_df['cond_code'] == cond) & (stats_df['selected'])].copy()
    if sub.empty:
        logging.warning('%s: no selected proteins; skipped.', COND_LABELS[cond])
        return
    sub = sub.sort_values('Chrom_FC', ascending=False).reset_index(drop=True)
    order = sub['Accession'].tolist()
    gene_of = dict(zip(sub['Accession'], sub['Gene symbol']))
    n_prot = len(order)

    base = f'{base_tag}_{COND_SHORT[cond]}'
    write_plotmatrix(outdir / f'{base}_plotmatrix.csv', order, gene_of, sub, [cond])

    chrom_row = sub['Chrom_FC'].fillna(0).values
    cyto_row  = sub['Cyto_FC'].fillna(0).values
    src_row   = sub['sig_source'].astype(int).values        # 1 chrom,2 cyto,3 both
    cand_row  = sub['transloc_candidate'].astype(int).values

    norm = TwoSlopeNorm(vmin=-args.vmax, vcenter=0, vmax=args.vmax)
    SRC_CMAP = ListedColormap(SRC_COLORS)
    SRC_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], SRC_CMAP.N)
    CAND_CMAP = ListedColormap(['#E8E8E8', '#2CA25F'])
    CAND_NORM = BoundaryNorm([-0.5, 0.5, 1.5], CAND_CMAP.N)

    FRAC_ROW_H = 0.55; BAND_H = 2 * FRAC_ROW_H
    ANNOT_ROW_H = 0.18; ANNOT_GAP = 0.10; STRIP_GAP = 0.18
    LEFT_IN = 2.70; RIGHT_IN = 2.55; TOP_PAD = 0.85; BOT_PAD = 1.25
    CBAR_W_IN = 2.8; CBAR_H_IN = 0.14; CBAR_BOT = 0.46

    hm_w = hm_width(n_prot)
    FIG_W = LEFT_IN + hm_w + RIGHT_IN
    FIG_H = TOP_PAD + 2 * ANNOT_ROW_H + ANNOT_GAP + STRIP_GAP + BAND_H + BOT_PAD

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    y_cursor = FIG_H - TOP_PAD

    # strip 1: significance source (chrom-only / cyto-only / both)
    y1 = y_cursor - ANNOT_ROW_H
    ax1 = fig.add_axes([LEFT_IN / FIG_W, y1 / FIG_H, hm_w / FIG_W, ANNOT_ROW_H / FIG_H])
    ax1.imshow(src_row.reshape(1, -1), aspect='auto', cmap=SRC_CMAP, norm=SRC_NORM,
               interpolation='nearest')
    ax1.set_xticks([]); ax1.set_yticks([])
    for sp in ax1.spines.values():
        sp.set_visible(False)
    fig.text((LEFT_IN - 0.12) / FIG_W, (y1 + ANNOT_ROW_H / 2) / FIG_H,
             'significance\nsource', ha='right', va='center',
             fontsize=7, color='#555555', style='italic')
    y_cursor = y1 - ANNOT_GAP

    # strip 2: translocation candidate (informational, not a filter)
    y2 = y_cursor - ANNOT_ROW_H
    ax2 = fig.add_axes([LEFT_IN / FIG_W, y2 / FIG_H, hm_w / FIG_W, ANNOT_ROW_H / FIG_H])
    ax2.imshow(cand_row.reshape(1, -1), aspect='auto', cmap=CAND_CMAP, norm=CAND_NORM,
               interpolation='nearest')
    ax2.set_xticks([]); ax2.set_yticks([])
    for sp in ax2.spines.values():
        sp.set_visible(False)
    fig.text((LEFT_IN - 0.12) / FIG_W, (y2 + ANNOT_ROW_H / 2) / FIG_H,
             'translocation\ncandidate', ha='right', va='center',
             fontsize=7, color='#555555', style='italic')
    y_cursor = y2 - STRIP_GAP

    y_band_bot = y_cursor - BAND_H
    im = draw_band(fig, FIG_W, FIG_H, LEFT_IN, y_band_bot, hm_w, FRAC_ROW_H,
                   chrom_row, cyto_row, norm)

    if n_prot <= 80:
        bottom_y = (y_band_bot - 0.04) / FIG_H
        for pi, acc in enumerate(order):
            xf = (LEFT_IN + (pi + 0.5) / n_prot * hm_w) / FIG_W
            fig.text(xf, bottom_y, gene_of.get(acc, acc), ha='right', va='top',
                     rotation=90, fontsize=5.5, color='#222222', style='italic')

    legend_patches = [mpatches.Patch(color=SRC_COLORS[1], label='sig: chromatin only'),
                      mpatches.Patch(color=SRC_COLORS[2], label='sig: cytosol only'),
                      mpatches.Patch(color=SRC_COLORS[3], label='sig: both'),
                      mpatches.Patch(color='#2CA25F',
                                     label='transloc. candidate (chrom↑cyto↓ / chrom↓cyto↑)')]
    fig.legend(handles=legend_patches, loc='lower right',
               bbox_to_anchor=(0.99, CBAR_BOT / FIG_H + 0.05),
               fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
               title='Annotation strips', title_fontsize=7)

    cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
    ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H, CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3]); cb.ax.tick_params(labelsize=7)
    ax_cb.set_xlabel('log2 fold change vs Naive', fontsize=8, color='#444444', labelpad=4)

    n_ch = int(sub['chrom_sig'].sum()); n_cy = int(sub['cyto_sig'].sum())
    n_bo = int(sub['both_sig'].sum()); n_cand = int(cand_row.sum())
    fig.text(0.5, (FIG_H - TOP_PAD * 0.30) / FIG_H,
             f'{COND_LABELS[cond]} vs Naive — significant proteins ({args.select_mode})  '
             f'(n = {n_prot};  chrom {n_ch}, cyto {n_cy}, both {n_bo};  candidates {n_cand})',
             ha='center', va='top', fontsize=11, fontweight='bold', color='#1A3A5C')
    fig.text(0.5, 0.008,
             f'Selection: |FC|≥{args.fc_thresh:g} & Corrected≥{args.corr_thresh:g} '
             f'in {("chromatin OR cytosol" if args.select_mode=="union" else args.select_mode)} '
             f'(up OR down) in {COND_LABELS[cond]}.  No ΔΔ filter.  '
             f'Chromatin filtered to {"/".join(args.keep_filter)}; cytosol unfiltered; AW-M-3 excluded.  '
             f'Sorted by Chrom_FC (up→down).',
             ha='center', va='bottom', fontsize=5.5, color='#888888')

    png = outdir / f'{base}.png'; pdf = outdir / f'{base}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logging.info('%s: wrote %s | %s (n=%d, candidates=%d)',
                 COND_LABELS[cond], png.name, pdf.name, n_prot, n_cand)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--corr-thresh', type=float, default=3.3)
    ap.add_argument('--fc-thresh', type=float, default=0.5)
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--keep-filter', nargs='+', default=['Keep', 'Review'])
    ap.add_argument('--select-mode', choices=['union', 'both', 'chrom', 'cyto'],
                    default='union',
                    help="Per-condition selection: union = sig in chromatin OR cytosol "
                         "(default); both = sig in BOTH (intersection); chrom/cyto = "
                         "that compartment only.")
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--file-tag', default='RatAlcoholProteome')
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    for noisy in ('matplotlib', 'fontTools', 'fontTools.subset', 'PIL'):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    logging.info('Loading %s', args.xlsx.name)
    ch = pd.read_excel(args.xlsx, sheet_name='Chromatin')
    cy = pd.read_excel(args.xlsx, sheet_name='Cytosol')
    for df in (ch, cy):
        df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

    n_before = len(ch)
    ch = ch[ch['Filter'].isin(args.keep_filter)].copy()
    logging.info('Chromatin filter %s: %d -> %d proteins', args.keep_filter, n_before, len(ch))

    ch = ch.drop_duplicates('Accession').set_index('Accession')
    cy = cy.drop_duplicates('Accession').set_index('Accession')
    common = [a for a in ch.index if a in cy.index]
    logging.info('Common accessions (chromatin-filtered ∩ cytosol): %d', len(common))
    ch = ch.loc[common]; cy = cy.loc[common]

    stats_df = build_stats(ch, cy, common)

    # ---- per-compartment significance (criterion A; direction NOT filtered) --
    stats_df['chrom_sig'] = (
        (stats_df['Chrom_FC'].abs() >= args.fc_thresh) &
        (stats_df['Chrom_Corrected'] >= args.corr_thresh)).fillna(False)
    stats_df['cyto_sig'] = (
        (stats_df['Cyto_FC'].abs() >= args.fc_thresh) &
        (stats_df['Cyto_Corrected'] >= args.corr_thresh)).fillna(False)
    stats_df['both_sig'] = stats_df['chrom_sig'] & stats_df['cyto_sig']
    stats_df['either_sig'] = stats_df['chrom_sig'] | stats_df['cyto_sig']

    # selection driving the figures (per condition), per --select-mode
    sel = {'union': stats_df['either_sig'], 'both': stats_df['both_sig'],
           'chrom': stats_df['chrom_sig'], 'cyto': stats_df['cyto_sig']}[args.select_mode]
    stats_df['selected'] = sel.fillna(False)

    # significance source label (only meaningful where selected): 0 none,1 chrom,2 cyto,3 both
    src = np.select(
        [stats_df['both_sig'], stats_df['chrom_sig'] & ~stats_df['cyto_sig'],
         stats_df['cyto_sig'] & ~stats_df['chrom_sig']], [3, 1, 2], default=0)
    stats_df['sig_source'] = src
    stats_df['sig_source_label'] = pd.Series(src).map(
        {0: 'none', 1: 'chrom-only', 2: 'cyto-only', 3: 'both'}).values

    # informational candidate flag: chrom up & cyto down, or chrom down & cyto up
    up_dn = (stats_df['Chrom_FC'] >= args.fc_thresh) & (stats_df['Cyto_FC'] <= -args.fc_thresh)
    dn_up = (stats_df['Chrom_FC'] <= -args.fc_thresh) & (stats_df['Cyto_FC'] >= args.fc_thresh)
    stats_df['transloc_candidate'] = (stats_df['selected'] & (up_dn | dn_up)).fillna(False)

    MODE_TAG = {'union': 'SigEitherCompartment', 'both': 'SigBothCompartments',
                'chrom': 'SigChromatin', 'cyto': 'SigCytosol'}
    base_tag = (f'{args.file_tag}_{MODE_TAG[args.select_mode]}_CytoChrom_Heatmap_'
                f'Corr{args.corr_thresh:g}_FC{args.fc_thresh:g}_{args.date}').replace('.', 'p')

    # round numeric columns for the saved full-stats CSV
    full = stats_df.copy()
    for c in ['Chrom_FC', 'Chrom_Corrected', 'Cyto_FC', 'Cyto_Corrected', 'interaction_ddelta']:
        full[c] = full[c].round(4)
    full_csv = args.outdir / f'{base_tag}_fullstats.csv'
    full.to_csv(full_csv, index=False)
    logging.info('Wrote full stats: %s (%d rows)', full_csv.name, len(full))

    # per-condition counts
    logging.info('Select mode: %s', args.select_mode)
    for cond in CONDITIONS:
        s = stats_df[(stats_df['cond_code'] == cond) & stats_df['selected']]
        n_ch = int(((stats_df['cond_code'] == cond) & stats_df['chrom_sig']).sum())
        n_cy = int(((stats_df['cond_code'] == cond) & stats_df['cyto_sig']).sum())
        n_bo = int(((stats_df['cond_code'] == cond) & stats_df['both_sig']).sum())
        cand = int(s['transloc_candidate'].sum())
        logging.info('%s: selected n=%d  (chrom_sig=%d, cyto_sig=%d, both=%d, candidates=%d)',
                     COND_LABELS[cond], len(s), n_ch, n_cy, n_bo, cand)
    n_union = stats_df.loc[stats_df['selected'], 'Accession'].nunique()
    logging.info('UNION over conditions (selected in ≥1): %d proteins', n_union)

    # ---- products ------------------------------------------------------------
    plot_union(stats_df, args, args.outdir, base_tag)
    for cond in CONDITIONS:
        plot_per_condition(stats_df, cond, args, args.outdir, base_tag)

    print('Done.')
    print('  full stats CSV:', full_csv.name)
    print('  union n:', n_union)


if __name__ == '__main__':
    main()
