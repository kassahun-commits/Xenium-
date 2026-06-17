#!/usr/bin/env python3
"""
Combined AW figure: pathway/GO enrichment (top) + true-translocation heatmap
(bottom), arranged in two direction columns so everything lines up.

  LEFT  column = Into Chromatin     (red)
  RIGHT column = Into Soluble Nuclear (blue)

  Panel A (top)    : Metascape GO/pathway enrichment bars (top-N terms),
                     extracted live from the two Metascape result PDFs.
  Panel B (bottom) : split heatmap of the top translocators per direction
                     (Chromatin & SN log2 fold-change vs Naive), computed
                     from the proteome workbook with the true-translocation
                     pipeline.

Standing rules (MEWS lab): workbook with AW-M-3 excluded; Chromatin filtered to
Filter in {Keep, Review}; SN uses ALL proteins; editable vector text
(pdf.fonttype=42); source CSVs alongside the figure; no hardcoded paths.

Usage:
  python3 make_combined_AW_metascape_heatmap_figure.py \
    --xlsx "../../EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx" \
    --up-pdf "/Users/naomi/Downloads/UP translcaotion.pdf" \
    --down-pdf "/Users/naomi/Downloads/down.pdf" \
    --outdir current --top-go 7 --top-per-dir 35 --date 2026-06-08
"""
import os
import argparse
import datetime as _dt
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize
from matplotlib.backends.backend_pdf import PdfPages
import fitz  # PyMuPDF

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

# ── sample columns (AW-M-3 already absent) ───────────────────────────────────
CH_N_COLS  = ['Chrom_N-F-1', 'Chrom_N-F-2', 'Chrom_N-F-3', 'Chrom_N-M-1', 'Chrom_N-M-2']
SN_N_COLS  = ['Nuc_N-F-1',   'Nuc_N-F-2',   'Nuc_N-F-3',   'Nuc_N-M-1',   'Nuc_N-M-2']
CH_AW_COLS = ['Chrom_AW-F-1', 'Chrom_AW-F-2', 'Chrom_AW-M-1', 'Chrom_AW-M-2']
SN_AW_COLS = ['Nuc_AW-F-1',  'Nuc_AW-F-2',  'Nuc_AW-M-1',  'Nuc_AW-M-2']

# ── heatmap diverging colormap (blue=SN, red=chromatin) ──────────────────────
HM_CMAP = LinearSegmentedColormap.from_list(
    'vp_diverging',
    ['#1A5FA0', '#2B7FD4', '#6AAEE0', '#AECFE8', '#DDEEF8',
     '#FFFFFF',
     '#FDD8E7', '#F5A0BC', '#EE5F8B', '#E8305A', '#C01E42'], N=512)

# ── shorten long GO term names (full names kept in source CSV) ───────────────
SHORTEN = {
    "alpha-amino acid metabolic process": "α-amino acid metabolism",
    "Metabolism of amino acids and derivatives": "Metabolism of amino acids & derivatives",
    "regulation of microtubule polymerization or depolymerization": "reg. of microtubule (de)polymerization",
    "aldehyde metabolic process": "aldehyde metabolism",
    "monocarboxylic acid metabolic process": "monocarboxylic acid metabolism",
    "positive regulation of cellular component biogenesis": "pos. reg. of cell component biogenesis",
    "regulation of microtubule cytoskeleton organization": "reg. of microtubule cytoskeleton org.",
    "purine ribonucleoside monophosphate metabolic process": "purine ribonucleoside-P metabolism",
    "regulation of plasma membrane bounded cell projection organization": "reg. of cell projection organization",
    "glutamine family amino acid metabolic process": "glutamine-family amino acid metabolism",
    "regulation of establishment of protein localization": "reg. of protein localization",
}


def shorten_name(name):
    return SHORTEN.get(name, name.replace("metabolic process", "metabolism"))


def go_cmap(direction):
    if direction == "up":
        return LinearSegmentedColormap.from_list(
            "ms_red", ["#FBD9E0", "#E8607F", "#C01E42", "#7A0B27"])
    return LinearSegmentedColormap.from_list(
        "ms_blue", ["#D6E5F3", "#5C97CE", "#1A5FA0", "#0C3460"])


# ── GO term extraction from Metascape PDF ────────────────────────────────────
def extract_terms(pdf_path):
    doc = fitz.open(pdf_path)
    pg = doc[0]
    words = pg.get_text("words")
    lines = defaultdict(list)
    for x0, y0, x1, y1, word, b, l, n in words:
        lines[(b, l)].append((x0, x1, y0, y1, word))

    term_rows, tick_words = [], []
    for ws in lines.values():
        ws = sorted(ws, key=lambda t: t[0])
        text = " ".join(t[4] for t in ws).strip()
        yc = sum((t[2] + t[3]) / 2 for t in ws) / len(ws)
        if ": " in text and any(text.startswith(p) for p in
                                ("GO:", "R-RNO", "rno", "WP", "hsa", "KEGG", "ko")):
            tid, name = text.split(": ", 1)
            term_rows.append((yc, tid.strip(), name.strip()))
        for x0, x1, y0, y1, word in ws:
            if word.strip().isdigit() and len(ws) == 1:
                tick_words.append(((x0 + x1) / 2, float(word)))

    xs = np.array([t[0] for t in tick_words])
    vs = np.array([t[1] for t in tick_words])
    a, b = np.polyfit(xs, vs, 1)

    bars = []
    for d in pg.get_drawings():
        fill = d.get("fill")
        if fill is None or min(fill) > 0.95:
            continue
        for item in d["items"]:
            if item[0] == "re":
                r = item[1]
                w, h = r.x1 - r.x0, r.y1 - r.y0
                if w > 2 and 4 < h < 16 and abs(r.x0 - 10.4) < 3:
                    bars.append(((r.y0 + r.y1) / 2, r.x1))
    out = []
    for yc, tid, name in term_rows:
        best = min(bars, key=lambda bb: abs(bb[0] - yc))
        if abs(best[0] - yc) > 8:
            continue
        out.append({"id": tid, "name": name, "value": float(a * best[1] + b)})
    out.sort(key=lambda d: d["value"], reverse=True)
    return out


# ── heatmap stats pipeline ───────────────────────────────────────────────────
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--up-pdf", required=True)
    ap.add_argument("--down-pdf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-go", type=int, default=7)
    ap.add_argument("--top-per-dir", type=int, default=35)
    ap.add_argument("--p-thresh", type=float, default=0.10)
    ap.add_argument("--fc-thresh", type=float, default=0.5)
    ap.add_argument("--vmax-hm", type=float, default=3.0)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- data ----
    up = extract_terms(args.up_pdf)[:args.top_go]
    down = extract_terms(args.down_pdf)[:args.top_go]
    go_vmax = max(t["value"] for t in up + down) * 1.02

    df = build_aw_stats(args.xlsx, args.p_thresh, args.fc_thresh)
    df_aw = df[df['TrueDir_AW'] != 'NS'].copy()
    n_ch_all = int((df_aw['TrueDir_AW'] == 'Into_Chromatin').sum())
    n_sn_all = int((df_aw['TrueDir_AW'] == 'Into_SN').sum())
    into_ch = (df_aw[df_aw['TrueDir_AW'] == 'Into_Chromatin']
               .sort_values('Interaction_AW', ascending=False)
               .head(args.top_per_dir).reset_index(drop=True))
    into_sn = (df_aw[df_aw['TrueDir_AW'] == 'Into_SN']
               .sort_values('Interaction_AW', ascending=True)
               .head(args.top_per_dir).reset_index(drop=True))
    print(f"GO terms: up={len(up)} down={len(down)}")
    print(f"AW true transloc: Into Chromatin={n_ch_all}, Into SN={n_sn_all}; "
          f"heatmap shows {len(into_ch)}/{len(into_sn)}")

    # ---- layout (inches) ----
    LEFT_PAD = 0.25
    GO_LABEL_W = 2.75
    GO_BAR_W = 2.35
    COL_W = GO_LABEL_W + GO_BAR_W       # 5.10
    COL_GAP = 0.60
    RIGHT_PAD = 0.25
    FIG_W = LEFT_PAD + 2 * COL_W + COL_GAP + RIGHT_PAD

    BOT_PAD = 1.15
    ROW_H = 0.16
    HM_H = args.top_per_dir * ROW_H
    MID_GAP = 1.00
    GO_H = 2.65
    TOP_PAD = 1.15
    FIG_H = BOT_PAD + HM_H + MID_GAP + GO_H + TOP_PAD

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    def ax_in(x, y, w, h):
        return fig.add_axes([x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H])

    # column x-origins
    xL0 = LEFT_PAD
    xR0 = LEFT_PAD + COL_W + COL_GAP
    # GO bar axes left edge (labels extend left of this)
    xL_bar = xL0 + GO_LABEL_W
    xR_bar = xR0 + GO_LABEL_W
    y_go = BOT_PAD + HM_H + MID_GAP
    # data-region centers (used to align heatmap under GO bars)
    cxL = xL_bar + GO_BAR_W / 2
    cxR = xR_bar + GO_BAR_W / 2

    # ---- Panel A: GO bar charts ----
    def draw_go(x_bar, terms, direction, vmax):
        ax = ax_in(x_bar, y_go, GO_BAR_W, GO_H)
        cmap = go_cmap(direction)
        norm = Normalize(vmin=0, vmax=vmax)
        n = len(terms)
        ypos = np.arange(n)[::-1]
        vals = [t["value"] for t in terms]
        names = [shorten_name(t["name"]) for t in terms]
        ax.barh(ypos, vals, height=0.72, color=[cmap(norm(v)) for v in vals], zorder=3)
        for y, v in zip(ypos, vals):
            ax.text(v + vmax * 0.02, y, f"{v:.1f}", va="center", ha="left",
                    fontsize=8, color="#222222", zorder=4)
        ax.set_yticks(ypos)
        ax.set_yticklabels(names, fontsize=9.5)
        ax.set_xlim(0, vmax * 1.14)
        ax.set_ylim(-0.7, n - 0.3)
        ax.set_xlabel(r"$-\log_{10}(P)$", fontsize=9)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", labelsize=8)
        ax.set_axisbelow(True)

    draw_go(xL_bar, up, "up", go_vmax)
    draw_go(xR_bar, down, "down", go_vmax)

    # ---- Panel B: heatmap blocks ----
    vmax = args.vmax_hm
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    FC_COL_W = 0.55
    BLOCK_DATA_W = 2 * FC_COL_W

    def draw_block(cx, sub):
        x0 = cx - FC_COL_W
        ax = ax_in(x0, BOT_PAD, BLOCK_DATA_W, HM_H)
        mat = np.column_stack([
            pd.to_numeric(sub['Chrom_FC_AW'], errors='coerce').fillna(0).values,
            pd.to_numeric(sub['SN_FC_AW'], errors='coerce').fillna(0).values])
        im = ax.imshow(mat, aspect='auto', cmap=HM_CMAP, norm=norm,
                       interpolation='nearest')
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub['Gene'].tolist(), fontsize=6.0)
        ax.tick_params(axis='y', length=0, pad=2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Chromatin', 'Soluble\nNuclear'], fontsize=6.5)
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', length=0, pad=3)
        ax.axvline(0.5, color='white', lw=1.2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        return im

    im = draw_block(cxL, into_ch)
    draw_block(cxR, into_sn)

    # ---- shared heatmap colorbar (bottom center) ----
    cbar_w = 2.6
    ax_cb = ax_in((FIG_W - cbar_w) / 2, 0.42, cbar_w, 0.12)
    cb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
    cb.set_ticks([-vmax, -2, -1, 0, 1, 2, vmax])
    cb.ax.tick_params(labelsize=6)
    ax_cb.set_xlabel('Log2 fold change vs Naive  (within fraction)',
                     fontsize=7, color='#444444', labelpad=3)

    # ---- headers, panel labels, title ----
    def fx(x):
        return x / FIG_W

    def fy(y):
        return y / FIG_H

    # overall title
    fig.text(0.5, fy(FIG_H - 0.30),
             "Acute Withdrawal — true Chromatin ↔ Soluble Nuclear translocation",
             ha="center", va="top", fontsize=14, fontweight="bold", color="#1A3A5C")

    # column headers (colored by direction)
    y_head = fy(y_go + GO_H + 0.34)
    fig.text(fx(cxL), y_head, "Into Chromatin", ha="center", va="bottom",
             fontsize=12.5, fontweight="bold", color="#C01E42")
    fig.text(fx(cxL), fy(y_go + GO_H + 0.16), f"({n_ch_all} proteins)",
             ha="center", va="bottom", fontsize=8, color="#C01E42")
    fig.text(fx(cxR), y_head, "Into Soluble Nuclear", ha="center", va="bottom",
             fontsize=12.5, fontweight="bold", color="#1A5FA0")
    fig.text(fx(cxR), fy(y_go + GO_H + 0.16), f"({n_sn_all} proteins)",
             ha="center", va="bottom", fontsize=8, color="#1A5FA0")

    # panel letters
    fig.text(fx(0.10), fy(y_go + GO_H + 0.10), "A", ha="left", va="bottom",
             fontsize=15, fontweight="bold")
    fig.text(fx(0.10), fy(BOT_PAD + HM_H + 0.30), "B", ha="left", va="bottom",
             fontsize=15, fontweight="bold")

    # panel sub-labels
    fig.text(fx(0.10), fy(y_go + GO_H - 0.02), "Pathway / GO enrichment",
             ha="left", va="top", fontsize=8.5, style="italic", color="#555555")
    fig.text(fx(0.10), fy(BOT_PAD + HM_H + 0.16),
             f"Top {args.top_per_dir} translocators per direction",
             ha="left", va="top", fontsize=8.5, style="italic", color="#555555")

    # footnote
    fig.text(0.5, fy(0.10),
             f"True translocation: BH FDR p_adj < {args.p_thresh}, |FC| ≥ {args.fc_thresh}, "
             f"opposite-direction movement  |  AW-M-3 excluded  |  "
             f"Chromatin: Keep+Review, SN: all proteins  |  "
             f"GO enrichment computed on full direction gene sets (Metascape)",
             ha="center", va="bottom", fontsize=5.2, color="#888888")

    # ---- save ----
    pthr = str(args.p_thresh).replace(".", "p")
    fcthr = str(args.fc_thresh).replace(".", "p")
    base = (f"{args.file_tag}_AW_CombinedFigure_GOenrich_SplitHeatmap"
            f"_top{args.top_go}GO_top{2*args.top_per_dir}HM_p{pthr}_FC{fcthr}_{args.date}")
    out_pdf = os.path.join(args.outdir, base + ".pdf")
    out_png = os.path.join(args.outdir, base + ".png")
    out_csv_hm = os.path.join(args.outdir, base + "_heatmap_source.csv")
    out_csv_go = os.path.join(args.outdir, base + "_GO_source.csv")

    pd.concat([into_ch, into_sn])[['Gene', 'Chrom_FC_AW', 'SN_FC_AW',
        'Interaction_AW', 'p_adj_AW', 'TrueDir_AW']].to_csv(out_csv_hm, index=False)
    go_rows = ([{"direction": "Into_Chromatin", "rank": i + 1, **t} for i, t in enumerate(up)] +
               [{"direction": "Into_SN", "rank": i + 1, **t} for i, t in enumerate(down)])
    pd.DataFrame(go_rows)[["direction", "rank", "id", "name", "value"]].to_csv(out_csv_go, index=False)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor="white")
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved:")
    for p in (out_pdf, out_png, out_csv_hm, out_csv_go):
        print(" ", p)


if __name__ == "__main__":
    main()
