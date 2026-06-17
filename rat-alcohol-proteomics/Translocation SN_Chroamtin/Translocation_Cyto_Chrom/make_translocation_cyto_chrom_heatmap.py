#!/usr/bin/env python3
"""
Cytosol -> Chromatin translocation heatmap, by condition (rat alcohol proteome).
===============================================================================
Two-compartment translocation model: CYTOSOL <-> CHROMATIN, replacing the older
SN <-> Chromatin model.  Soluble-nuclear is NOT used here.

A protein is called a TRUE cytosol->chromatin translocation in a condition
(Intoxication / Acute Withdrawal / Protracted Abstinence, each vs Naive) when
ALL THREE hold:

  1. Up & significant on chromatin :  Chrom_FC >= +FC_THRESH  AND
                                      Chrom_Corrected >= CORR_THRESH
  2. Down in cytosol              :  Cyto_FC  <= -FC_THRESH
  3. Genuine redistribution       :  delta-delta interaction test significant
                                     (p_adj < P_THRESH) AND interaction > 0
                                     (i.e. the chromatin gain is significantly
                                     larger than the cytosol change)

delta-delta interaction (recomputed from per-sample log2 abundances):
    delta_chrom = chrom_cond - mean(chrom_naive)   (per replicate)
    delta_cyto  = cyto_cond  - mean(cyto_naive)     (per replicate)
    interaction = mean(delta_chrom) - mean(delta_cyto)
    p           = Welch t-test(delta_chrom, delta_cyto), BH-corrected per condition

Fold change (display + criterion #1/#2) is log2 vs Naive, taken from the sheet's
precomputed 'Fold change' columns (verified == mean(cond)-mean(naive) of the
AW-M-3-excluded per-sample columns).  'Corrected' is the sheet's significance
score (tracks -log10 p; Corrected~=4.0 at p~=0.05, so 3.3 ~ p<=0.07).

Filtering: CHROMATIN restricted to Filter in {Keep, Review}.  Cytosol is NOT
filtered.  AW-M-3 is already excluded from the input workbook's fraction sheets.

MEWS rules: no hardcoded paths (CLI), editable vector text (pdf.fonttype=42),
source-data CSV (full per-protein stats + the plotted matrix) alongside figures.
Outputs: PNG (300 dpi), editable PDF, two CSVs (full stats + plotted matrix).
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
# sheet's precomputed log2 fold-change + 'Corrected' significance, per condition
SHEET_FC   = {'I': 'Fold change',   'AW': 'Fold change.1', 'PA': 'Fold change.2'}
SHEET_CORR = {'I': 'Corrected',     'AW': 'Corrected.1',   'PA': 'Corrected.2'}

CONDITIONS  = ['I', 'AW', 'PA']
COND_LABELS = {'I': 'Intoxication', 'AW': 'Acute\nWithdrawal', 'PA': 'Protracted\nAbstinence'}

# diverging colormap (blue = down vs naive, red = up) -- matches prior figures
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True, type=Path,
                    help='Path to the "Excluding AWM3" proteome workbook.')
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--corr-thresh', type=float, default=3.3,
                    help="Chromatin 'Corrected' significance cutoff (default 3.3).")
    ap.add_argument('--fc-thresh', type=float, default=0.5,
                    help='|log2 fold change| cutoff (default 0.5).')
    ap.add_argument('--p-thresh', type=float, default=0.10,
                    help='BH-adjusted p cutoff for the delta-delta test (default 0.10).')
    ap.add_argument('--vmax', type=float, default=3.0, help='Colorbar |max| (default 3).')
    ap.add_argument('--keep-filter', nargs='+', default=['Keep', 'Review'],
                    help="Chromatin 'Filter' values to retain (default: Keep Review).")
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--file-tag', default='RatAlcoholProteome')
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')

    # ---- load ----------------------------------------------------------------
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
    common = [a for a in ch.index if a in cy.index]          # chromatin universe ∩ cytosol
    logging.info('Common accessions (chromatin-filtered ∩ cytosol): %d', len(common))
    ch = ch.loc[common]
    cy = cy.loc[common]

    # ---- per-condition stats -------------------------------------------------
    ch_naive = vals(ch, CH_SAMPLES['N'])
    cy_naive = vals(cy, CY_SAMPLES['N'])
    ch_naive_mean = np.nanmean(ch_naive, axis=1)
    cy_naive_mean = np.nanmean(cy_naive, axis=1)

    rows = []
    pmat = {}
    for cond in CONDITIONS:
        ch_c = vals(ch, CH_SAMPLES[cond])
        cy_c = vals(cy, CY_SAMPLES[cond])
        d_ch = ch_c - ch_naive_mean[:, None]                 # per-replicate deltas
        d_cy = cy_c - cy_naive_mean[:, None]
        interaction = np.nanmean(d_ch, axis=1) - np.nanmean(d_cy, axis=1)

        pvec = np.full(len(common), np.nan)
        for i in range(len(common)):
            a = d_ch[i][np.isfinite(d_ch[i])]
            b = d_cy[i][np.isfinite(d_cy[i])]
            if a.size >= 2 and b.size >= 2:
                pvec[i] = stats.ttest_ind(a, b, equal_var=False).pvalue
        padj = bh_correction(pvec)
        pmat[cond] = (interaction, pvec, padj)

        chrom_fc = pd.to_numeric(ch[SHEET_FC[cond]], errors='coerce').values
        cyto_fc  = pd.to_numeric(cy[SHEET_FC[cond]], errors='coerce').values
        chrom_corr = pd.to_numeric(ch[SHEET_CORR[cond]], errors='coerce').values

        for i, acc in enumerate(common):
            rows.append({
                'Accession': acc,
                'Gene symbol': ch['Gene symbol'].iloc[i],
                'condition': COND_LABELS[cond].replace('\n', ' '),
                'cond_code': cond,
                'Chrom_FC': round(float(chrom_fc[i]), 4) if np.isfinite(chrom_fc[i]) else np.nan,
                'Chrom_Corrected': round(float(chrom_corr[i]), 4) if np.isfinite(chrom_corr[i]) else np.nan,
                'Cyto_FC': round(float(cyto_fc[i]), 4) if np.isfinite(cyto_fc[i]) else np.nan,
                'interaction_ddelta': round(float(interaction[i]), 4) if np.isfinite(interaction[i]) else np.nan,
                'welch_p': pvec[i],
                'welch_p_adj': padj[i],
            })

    stats_df = pd.DataFrame(rows)

    # ---- true-translocation rule (cytosol -> chromatin) ---------------------
    def is_hit(r):
        return bool(
            np.isfinite(r['Chrom_FC']) and r['Chrom_FC'] >= args.fc_thresh and
            np.isfinite(r['Chrom_Corrected']) and r['Chrom_Corrected'] >= args.corr_thresh and
            np.isfinite(r['Cyto_FC']) and r['Cyto_FC'] <= -args.fc_thresh and
            np.isfinite(r['welch_p_adj']) and r['welch_p_adj'] < args.p_thresh and
            np.isfinite(r['interaction_ddelta']) and r['interaction_ddelta'] > 0
        )
    stats_df['translocation_hit'] = stats_df.apply(is_hit, axis=1)

    hits_by_cond = stats_df[stats_df['translocation_hit']].groupby('cond_code')['Accession'].apply(set)
    for cond in CONDITIONS:
        n = len(hits_by_cond.get(cond, set()))
        logging.info('%s: cytosol->chromatin hits = %d', COND_LABELS[cond].replace(chr(10), ' '), n)

    hit_acc = set(stats_df.loc[stats_df['translocation_hit'], 'Accession'])
    logging.info('Union (hit in >=1 condition): %d proteins', len(hit_acc))

    # ---- build plotted matrix (union proteins) ------------------------------
    union = sorted(hit_acc)
    base = (f'{args.file_tag}_Translocation_CytoToChrom_byCondition_Heatmap_'
            f'Corr{args.corr_thresh:g}_FC{args.fc_thresh:g}_p{args.p_thresh:g}_{args.date}'
            ).replace('.', 'p')

    # always write the full stats table (every tested protein)
    stats_csv = args.outdir / f'{base}_fullstats.csv'
    stats_df.to_csv(stats_csv, index=False)

    if not union:
        logging.warning('No proteins met the translocation rule; figure not drawn.')
        print('Wrote (no hits):', stats_csv.name)
        return

    sub = stats_df[stats_df['Accession'].isin(union)].copy()
    gene_of = dict(zip(sub['Accession'], sub['Gene symbol']))

    def mat(field):
        m = sub.pivot(index='Accession', columns='cond_code', values=field).reindex(union)
        return m[CONDITIONS]

    chrom_mat = mat('Chrom_FC')
    cyto_mat  = mat('Cyto_FC')
    hit_mat   = mat('translocation_hit').fillna(False).astype(bool)

    # sort proteins by total chromatin gain (strongest translocators first)
    order = chrom_mat.fillna(0).sum(axis=1).sort_values(ascending=False).index.tolist()
    chrom_mat = chrom_mat.reindex(order); cyto_mat = cyto_mat.reindex(order)
    hit_mat = hit_mat.reindex(order)
    n_prot = len(order)

    # plotted-matrix CSV (exact values behind the figure)
    plot_csv = args.outdir / f'{base}_plotmatrix.csv'
    out_rows = []
    for acc in order:
        for cond in CONDITIONS:
            out_rows.append({'Accession': acc, 'Gene symbol': gene_of.get(acc, ''),
                             'condition': COND_LABELS[cond].replace('\n', ' '),
                             'Chrom_FC': chrom_mat.loc[acc, cond],
                             'Cyto_FC': cyto_mat.loc[acc, cond],
                             'translocation_hit': bool(hit_mat.loc[acc, cond])})
    pd.DataFrame(out_rows).to_csv(plot_csv, index=False)

    # ---- figure --------------------------------------------------------------
    norm = TwoSlopeNorm(vmin=-args.vmax, vcenter=0, vmax=args.vmax)
    ANNOT_CMAP = ListedColormap(['#E8E8E8', '#E8305A'])      # 0 = NS, 1 = hit
    ANNOT_NORM = BoundaryNorm([-0.5, 0.5, 1.5], ANNOT_CMAP.N)

    FRAC_ROW_H = 0.50; BAND_H = 2 * FRAC_ROW_H
    ANNOT_ROW_H = 0.16; ANNOT_H = 3 * ANNOT_ROW_H
    ANNOT_GAP = 0.20; SECTION_GAP = 0.30; N_COND = 3
    LEFT_IN = 2.70; RIGHT_IN = 2.90; TOP_PAD = 0.60; BOT_PAD = 1.25
    CBAR_W_IN = 2.8; CBAR_H_IN = 0.14; CBAR_BOT = 0.46

    hm_w = max(n_prot * 0.10, 3.0)
    FIG_W = LEFT_IN + hm_w + RIGHT_IN
    FIG_H = (TOP_PAD + ANNOT_H + ANNOT_GAP
             + N_COND * BAND_H + (N_COND - 1) * SECTION_GAP + BOT_PAD)

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
    y_cursor = FIG_H - TOP_PAD

    # annotation strip: hit / NS per condition
    annot = np.vstack([hit_mat[c].astype(int).values for c in CONDITIONS])
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
        fig.text((LEFT_IN + hm_w + 0.10) / FIG_W, yf, COND_LABELS[cond].replace('\n', ' '),
                 ha='left', va='center', fontsize=6, color='#333333')
    fig.text((LEFT_IN - 0.12) / FIG_W, (y_annot_bot + ANNOT_H / 2) / FIG_H,
             'Cytosol→Chromatin\ntranslocation', ha='right', va='center',
             fontsize=7, color='#555555', style='italic')
    y_cursor = y_annot_bot - ANNOT_GAP

    COMP_BG = {'Chromatin': '#FBF9F0', 'Cytosol': '#F0F6FB'}
    im = None
    for ci, cond in enumerate(CONDITIONS):
        y_band_bot = y_cursor - BAND_H

        ch_row = chrom_mat[cond].fillna(0).values.reshape(1, -1)
        ax_ch = fig.add_axes([LEFT_IN / FIG_W, (y_band_bot + FRAC_ROW_H) / FIG_H,
                              hm_w / FIG_W, FRAC_ROW_H / FIG_H], facecolor=COMP_BG['Chromatin'])
        im = ax_ch.imshow(ch_row, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
        ax_ch.set_xticks([]); ax_ch.set_yticks([])
        for sp in ax_ch.spines.values():
            sp.set_visible(False)
        fig.text((LEFT_IN + hm_w + 0.10) / FIG_W,
                 (y_band_bot + FRAC_ROW_H + FRAC_ROW_H / 2) / FIG_H, 'Chromatin',
                 ha='left', va='center', fontsize=7.5, color='#111111')

        cy_row = cyto_mat[cond].fillna(0).values.reshape(1, -1)
        ax_cy = fig.add_axes([LEFT_IN / FIG_W, y_band_bot / FIG_H,
                              hm_w / FIG_W, FRAC_ROW_H / FIG_H], facecolor=COMP_BG['Cytosol'])
        ax_cy.imshow(cy_row, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
        ax_cy.set_xticks([]); ax_cy.set_yticks([])
        for sp in ax_cy.spines.values():
            sp.set_visible(False)
        fig.text((LEFT_IN + hm_w + 0.10) / FIG_W, (y_band_bot + FRAC_ROW_H / 2) / FIG_H,
                 'Cytosol', ha='left', va='center', fontsize=7.5, color='#111111')

        fig.text((LEFT_IN - 0.12) / FIG_W, (y_band_bot + BAND_H / 2) / FIG_H,
                 COND_LABELS[cond].replace('\n', ' '), ha='right', va='center',
                 fontsize=10, fontweight='bold', color='#111111')
        sep_y = (y_band_bot + FRAC_ROW_H) / FIG_H
        fig.add_artist(plt.Line2D([LEFT_IN / FIG_W, (LEFT_IN + hm_w) / FIG_W], [sep_y, sep_y],
                                  transform=fig.transFigure, color='#CCCCCC', linewidth=0.8, zorder=5))
        y_cursor = y_band_bot - SECTION_GAP

    # gene labels along the bottom (vertical) when not too many
    if n_prot <= 80:
        bottom_y = (y_cursor + SECTION_GAP - 0.02) / FIG_H
        for pi, acc in enumerate(order):
            xf = (LEFT_IN + (pi + 0.5) / n_prot * hm_w) / FIG_W
            fig.text(xf, bottom_y, gene_of.get(acc, acc), ha='right', va='top',
                     rotation=90, fontsize=5.5, color='#222222', style='italic')

    legend_patches = [
        mpatches.Patch(color='#E8305A', label='Cytosol→Chromatin hit'),
        mpatches.Patch(color='#E8E8E8', label='not a hit', linewidth=0.5, edgecolor='#AAAAAA'),
    ]
    fig.legend(handles=legend_patches, loc='lower right',
               bbox_to_anchor=(0.99, CBAR_BOT / FIG_H + 0.05),
               fontsize=7, framealpha=0.9, edgecolor='#CCCCCC',
               title='Annotation strip', title_fontsize=7)

    cbar_left = (FIG_W - CBAR_W_IN) / 2 / FIG_W
    ax_cb = fig.add_axes([cbar_left, CBAR_BOT / FIG_H, CBAR_W_IN / FIG_W, CBAR_H_IN / FIG_H])
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3]); cb.ax.tick_params(labelsize=7)
    ax_cb.set_xlabel('log2 fold change vs Naive', fontsize=8, color='#444444', labelpad=4)

    fig.text(0.5, (FIG_H - TOP_PAD * 0.35) / FIG_H,
             f'Cytosol → Chromatin translocation by condition  (n = {n_prot})',
             ha='center', va='top', fontsize=11, fontweight='bold', color='#1A3A5C')
    fig.text(0.5, 0.008,
             f'Hit = Chrom_FC≥+{args.fc_thresh:g} & Chrom Corrected≥{args.corr_thresh:g}  +  '
             f'Cyto_FC≤−{args.fc_thresh:g}  +  ΔΔ interaction p_adj<{args.p_thresh:g} (into chromatin).  '
             f'Chromatin filtered to {"/".join(args.keep_filter)}; cytosol unfiltered; AW-M-3 excluded.  '
             f'Union ≥1 condition.',
             ha='center', va='bottom', fontsize=5.5, color='#888888')

    png = args.outdir / f'{base}.png'
    pdf = args.outdir / f'{base}.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', facecolor='white')
    plt.close(fig)

    logging.info('Wrote %s | %s', png.name, pdf.name)
    logging.info('Wrote %s | %s', stats_csv.name, plot_csv.name)
    print('Done. Union hits:', n_prot, '| genes:',
          ', '.join(gene_of.get(a, a) for a in order))


if __name__ == '__main__':
    main()
