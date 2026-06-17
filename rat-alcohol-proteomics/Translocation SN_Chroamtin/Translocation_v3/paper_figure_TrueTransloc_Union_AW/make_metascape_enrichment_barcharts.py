#!/usr/bin/env python3
"""
Recreate Metascape enrichment bar charts (cleaner / paper-ready) for the
AW true-translocation gene sets.

  UP   = proteins translocating INTO CHROMATIN during Acute Withdrawal  (red)
  DOWN = proteins translocating INTO SOLUBLE NUCLEAR (SN) during AW      (blue)

The script does NOT hardcode the enrichment values: it re-extracts term names
and bar lengths directly from the two Metascape result PDFs at runtime using
PyMuPDF (fitz), calibrating the x-axis from the numeric tick labels so that
each bar's right edge maps to -log10(P).

Outputs (per MEWS lab rules): editable vector text (pdf.fonttype=42), a source
CSV alongside each figure, and a combined side-by-side panel figure.

Example:
  python3 make_metascape_enrichment_barcharts.py \
    --up-pdf "/Users/naomi/Downloads/UP translcaotion.pdf" \
    --down-pdf "/Users/naomi/Downloads/down.pdf" \
    --outdir current --top-n 15 --date 2026-06-08
"""
import argparse
import os
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

import fitz  # PyMuPDF

# ---- editable vector text ----
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"


def extract_terms(pdf_path):
    """Return list of dicts {id, name, value, color_frac} sorted most->least
    significant, extracted from a Metascape horizontal-bar PDF."""
    doc = fitz.open(pdf_path)
    pg = doc[0]
    words = pg.get_text("words")  # x0,y0,x1,y1,word,block,line,wordno

    # 1) group words into text lines by (block,line)
    lines = defaultdict(list)
    for x0, y0, x1, y1, word, b, l, n in words:
        lines[(b, l)].append((x0, x1, y0, y1, word))

    term_rows = []   # (yc, id, name)
    tick_words = []  # (xc, value) for single-digit numeric ticks
    for key, ws in lines.items():
        ws = sorted(ws, key=lambda t: t[0])
        text = " ".join(t[4] for t in ws).strip()
        yc = sum((t[2] + t[3]) / 2 for t in ws) / len(ws)
        # term lines look like  "GO:0034504: protein localization to nucleus"
        if ":" in text and any(text.startswith(p) for p in ("GO:", "R-RNO", "rno", "WP", "hsa", "KEGG", "ko")):
            # split id from name on the first ": " after the accession
            if ": " in text:
                tid, name = text.split(": ", 1)
            else:
                tid, name = text, text
            term_rows.append((yc, tid.strip(), name.strip()))
        # numeric tick labels: a single token that is an integer
        for x0, x1, y0, y1, word in ws:
            if word.strip().isdigit() and len(ws) == 1:
                tick_words.append(((x0 + x1) / 2, float(word)))

    # 2) calibrate x -> value from tick labels (linear fit)
    if len(tick_words) < 2:
        raise RuntimeError(f"Could not find axis tick labels in {pdf_path}")
    xs = np.array([t[0] for t in tick_words])
    vs = np.array([t[1] for t in tick_words])
    a, b = np.polyfit(xs, vs, 1)  # value = a*x + b

    # 3) collect colored bar rectangles
    bars = []  # (yc, x1, fill)
    for d in pg.get_drawings():
        fill = d.get("fill")
        if fill is None:
            continue
        # skip near-white (page / plot background)
        if min(fill) > 0.95:
            continue
        for item in d["items"]:
            if item[0] == "re":
                r = item[1]
                w, h = r.x1 - r.x0, r.y1 - r.y0
                if w > 2 and 4 < h < 16 and abs(r.x0 - 10.4) < 3:
                    bars.append(((r.y0 + r.y1) / 2, r.x1, tuple(fill)))

    # 4) pair each term row to the bar at the same y
    out = []
    for yc, tid, name in term_rows:
        best = min(bars, key=lambda bb: abs(bb[0] - yc))
        if abs(best[0] - yc) > 8:
            continue
        value = a * best[1] + b
        out.append({"id": tid, "name": name, "value": float(value)})

    out.sort(key=lambda d: d["value"], reverse=True)
    return out


# --- shorten long term names for display (full name kept in source CSV) ---
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
    if name in SHORTEN:
        return SHORTEN[name]
    return name.replace("metabolic process", "metabolism")


def make_cmap(direction):
    if direction == "up":   # into chromatin -> reds
        return LinearSegmentedColormap.from_list(
            "ms_red", ["#FBD9E0", "#E8607F", "#C01E42", "#7A0B27"])
    else:                   # into SN -> blues
        return LinearSegmentedColormap.from_list(
            "ms_blue", ["#D6E5F3", "#5C97CE", "#1A5FA0", "#0C3460"])


def draw_chart(ax, terms, direction, vmin, vmax, title):
    cmap = make_cmap(direction)
    norm = Normalize(vmin=vmin, vmax=vmax)
    n = len(terms)
    ypos = np.arange(n)[::-1]  # most significant at top
    vals = [t["value"] for t in terms]
    names = [shorten_name(t["name"]) for t in terms]
    colors = [cmap(norm(v)) for v in vals]

    ax.barh(ypos, vals, height=0.74, color=colors, zorder=3)

    # value annotations at bar ends
    for y, v in zip(ypos, vals):
        ax.text(v + vmax * 0.015, y, f"{v:.1f}", va="center", ha="left",
                fontsize=19, color="#222222", zorder=4)

    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=24)
    ax.set_xlim(0, vmax * 1.12)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(r"$-\log_{10}(P)$", fontsize=22)
    ax.set_title(title, fontsize=27, fontweight="bold", pad=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=18)
    ax.set_axisbelow(True)
    return cmap, norm


def save_single(terms, direction, title, vmin, vmax, outpath_base):
    n = len(terms)
    fig_h = 0.85 * n + 1.8
    fig = plt.figure(figsize=(14, fig_h))
    ax = fig.add_axes([0.46, 0.10 if n > 8 else 0.16, 0.50, 0.80])
    draw_chart(ax, terms, direction, vmin, vmax, title)

    fig.savefig(outpath_base + ".pdf", bbox_inches="tight")
    fig.savefig(outpath_base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_combined(up, down, vmin, vmax, outpath_base):
    n = max(len(up), len(down))
    fig_h = 0.85 * n + 2.0
    fig = plt.figure(figsize=(26, fig_h))
    axL = fig.add_axes([0.235, 0.12, 0.235, 0.80])
    axR = fig.add_axes([0.745, 0.12, 0.235, 0.80])
    draw_chart(axL, up, "up", vmin, vmax, "Into Chromatin (AW)")
    draw_chart(axR, down, "down", vmin, vmax, "Into Soluble Nuclear (AW)")
    fig.savefig(outpath_base + ".pdf", bbox_inches="tight")
    fig.savefig(outpath_base + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_csv(terms, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "term_id", "term_name", "neg_log10_P"])
        for i, t in enumerate(terms, 1):
            w.writerow([i, t["id"], t["name"], f"{t['value']:.4f}"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--up-pdf", required=True)
    ap.add_argument("--down-pdf", required=True)
    ap.add_argument("--outdir", default="current")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--date", default="2026-06-08")
    ap.add_argument("--file-tag", default="")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    up_all = extract_terms(args.up_pdf)
    down_all = extract_terms(args.down_pdf)
    print(f"UP   terms extracted: {len(up_all)}")
    print(f"DOWN terms extracted: {len(down_all)}")

    n = args.top_n
    up = up_all[:n]
    down = down_all[:n]

    # shared color/axis scale so the two charts are directly comparable
    allv = [t["value"] for t in up + down]
    vmax = max(allv) * 1.02
    vmin = 0.0

    tag = f"_{args.file_tag}" if args.file_tag else ""
    base = os.path.join(args.outdir,
                        f"RatAlcoholProteome_AW_MetascapeEnrichment")

    up_base = f"{base}_IntoChromatin_top{n}{tag}_{args.date}"
    dn_base = f"{base}_IntoSN_top{n}{tag}_{args.date}"
    cb_base = f"{base}_Combined_top{n}{tag}_{args.date}"

    save_single(up, "up", f"Into Chromatin (AW) — top {n} terms",
                vmin, vmax, up_base)
    save_single(down, "down", f"Into Soluble Nuclear (AW) — top {n} terms",
                vmin, vmax, dn_base)
    save_combined(up, down, vmin, vmax, cb_base)

    write_csv(up, up_base + "_source.csv")
    write_csv(down, dn_base + "_source.csv")

    print("\nWrote:")
    for p in (up_base, dn_base, cb_base):
        print(" ", p + ".pdf")
    print("  source CSVs alongside each single chart.")


if __name__ == "__main__":
    main()
