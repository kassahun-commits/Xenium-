#!/usr/bin/env python3
"""
Cross-condition Venn diagrams of TRUE translocators (Chromatin <-> Soluble
Nuclear), one per direction, styled like the reference proportional Venn.

  LEFT  : Into Chromatin   — overlap of true into-chromatin translocators
          across AW / Intoxication / PA
  RIGHT : Into Soluble Nuclear — overlap of true into-SN translocators

True-translocation criteria (condition X vs Naive), identical to Panels A/B:
  Into Chromatin : Sig(p_adj<P) AND interaction>0 AND SN_FC   < -FC
  Into SN        : Sig(p_adj<P) AND interaction<0 AND SN_FC > 0 AND Chrom_FC < -FC
where interaction = mean(delta_chromatin) - mean(delta_SN), deltas are
condition replicates minus the naive mean (within each fraction), p from a
Welch t-test on the two delta vectors, BH-FDR corrected across common proteins
(per condition).

Standing rules (MEWS lab): workbook with AW-M-3 excluded; Chromatin filtered to
Filter in {Keep, Review}; SN uses ALL proteins; editable vector text
(pdf.fonttype=42); source CSVs alongside the figure; no hardcoded paths.

Usage:
  python3 make_crosscondition_translocation_venns.py \
    --xlsx "/.../EDIT  Excluding AWM3 ... copy.xlsx" \
    --outdir current --p-thresh 0.10 --fc-thresh 0.5 --date 2026-06-08
"""
import os
import argparse
import datetime as _dt

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn3, venn3_circles

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

# reference circle colors
COL_AW, COL_INTOX, COL_PA = "#5FAD8E", "#D4845A", "#9B8EC4"

# condition -> replicate suffixes (AW-M-3 absent => AW has 4 reps)
REPS = {
    "N":  ["F-1", "F-2", "F-3", "M-1", "M-2"],
    "I":  ["F-1", "F-2", "F-3", "M-1", "M-2"],
    "AW": ["F-1", "F-2", "M-1", "M-2"],
    "PA": ["F-1", "F-2", "F-3", "M-1", "M-2"],
}


def cols(prefix, cond):
    return [f"{prefix}{cond}-{r}" for r in REPS[cond]]


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


def get_vals(df, gene, c):
    if gene in df.index:
        return pd.to_numeric(df.loc[gene, c], errors="coerce").values
    return np.full(len(c), np.nan)


def mean_fc(df, gene, cond_c, naive_c):
    return float(np.nanmean(get_vals(df, gene, cond_c)) -
                 np.nanmean(get_vals(df, gene, naive_c)))


def build_cond_stats(ch_df, sn_df, common, cond, p_thresh, fc_thresh):
    ch_naive_cols, sn_naive_cols = cols("Chrom_", "N"), cols("Nuc_", "N")
    ch_cond_cols, sn_cond_cols = cols("Chrom_", cond), cols("Nuc_", cond)
    rows = []
    for g in common:
        d_ch = get_vals(ch_df, g, ch_cond_cols) - np.nanmean(get_vals(ch_df, g, ch_naive_cols))
        d_sn = get_vals(sn_df, g, sn_cond_cols) - np.nanmean(get_vals(sn_df, g, sn_naive_cols))
        interaction = float(np.nanmean(d_ch) - np.nanmean(d_sn))
        a = d_ch[np.isfinite(d_ch)]
        b = d_sn[np.isfinite(d_sn)]
        p = stats.ttest_ind(a, b, equal_var=False)[1] if len(a) >= 2 and len(b) >= 2 else np.nan
        rows.append({"Gene": g, "Interaction": interaction, "p": p,
                     "Chrom_FC": mean_fc(ch_df, g, ch_cond_cols, ch_naive_cols),
                     "SN_FC": mean_fc(sn_df, g, sn_cond_cols, sn_naive_cols)})
    df = pd.DataFrame(rows)
    valid = df["p"].notna()
    df["p_adj"] = np.nan
    df.loc[valid, "p_adj"] = bh_correction(df.loc[valid, "p"].values)
    sig = df["p_adj"] < p_thresh
    inter, snfc, chfc = df["Interaction"], df["SN_FC"], df["Chrom_FC"]
    true_ch = sig & (inter > 0) & (snfc < -fc_thresh)
    true_sn = sig & (inter < 0) & (snfc > 0) & (chfc < -fc_thresh)
    df["TrueDir"] = np.where(true_ch, "Into_Chromatin",
                    np.where(true_sn, "Into_SN", "NS"))
    return df


def draw_venn(ax, sets, totals, title):
    aw, intox, pa = sets
    v = venn3([aw, intox, pa], set_labels=(
        f"AW\n(n={totals[0]})", f"Intox.\n(n={totals[1]})", f"PA\n(n={totals[2]})"),
        set_colors=(COL_AW, COL_INTOX, COL_PA), alpha=0.85, ax=ax)
    c = venn3_circles([aw, intox, pa], linewidth=1.0, color="#5A5A5A", ax=ax)
    # style subset count labels
    for sid in ("100", "010", "001", "110", "101", "011", "111"):
        lbl = v.get_label_by_id(sid)
        if lbl:
            lbl.set_fontsize(13)
            lbl.set_color("#222222")
    # style set labels
    for sid in ("A", "B", "C"):
        lbl = v.get_label_by_id(sid)
        if lbl:
            lbl.set_fontsize(12.5)
            lbl.set_fontweight("bold")
    ax.set_title(title, fontsize=13, fontweight="bold", color="#1A3A5C", pad=10)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--p-thresh", type=float, default=0.10)
    ap.add_argument("--fc-thresh", type=float, default=0.5)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ch_raw = pd.read_excel(args.xlsx, sheet_name="Chromatin")
    ch_df = ch_raw[ch_raw["Filter"].isin(["Keep", "Review"])].copy()
    sn_df = pd.read_excel(args.xlsx, sheet_name="Soluble nuclear").copy()
    for d in (ch_df, sn_df):
        d["Gene symbol"] = d["Gene symbol"].astype(str).str.strip()
    ch_df = ch_df.drop_duplicates("Gene symbol").set_index("Gene symbol")
    sn_df = sn_df.drop_duplicates("Gene symbol").set_index("Gene symbol")
    common = sorted(set(ch_df.index) & set(sn_df.index))
    print(f"Common proteins: {len(common)}")

    cond_dir = {}      # cond -> {gene: TrueDir}
    for cond in ("I", "AW", "PA"):
        df = build_cond_stats(ch_df, sn_df, common, cond, args.p_thresh, args.fc_thresh)
        cond_dir[cond] = df.set_index("Gene")["TrueDir"]
        nch = int((df["TrueDir"] == "Into_Chromatin").sum())
        nsn = int((df["TrueDir"] == "Into_SN").sum())
        print(f"  {cond}: Into Chromatin={nch}, Into SN={nsn}")

    def sets_for(direction):
        return {cond: set(s.index[s == direction]) for cond, s in cond_dir.items()}

    ch_sets = sets_for("Into_Chromatin")
    sn_sets = sets_for("Into_SN")

    fig = plt.figure(figsize=(13, 6.4), facecolor="white")
    axL = fig.add_axes([0.02, 0.07, 0.46, 0.80])
    axR = fig.add_axes([0.52, 0.07, 0.46, 0.80])

    draw_venn(axL, (ch_sets["AW"], ch_sets["I"], ch_sets["PA"]),
              (len(ch_sets["AW"]), len(ch_sets["I"]), len(ch_sets["PA"])),
              "Into Chromatin — true translocators\noverlap across conditions")
    draw_venn(axR, (sn_sets["AW"], sn_sets["I"], sn_sets["PA"]),
              (len(sn_sets["AW"]), len(sn_sets["I"]), len(sn_sets["PA"])),
              "Into Soluble Nuclear — true translocators\noverlap across conditions")

    fig.text(0.5, 0.965, "True protein translocation — cross-condition overlap (Chromatin vs Soluble Nuclear)",
             ha="center", va="top", fontsize=14, fontweight="bold", color="#1A3A5C")
    fig.text(0.5, 0.02,
             f"True translocation: BH FDR p_adj < {args.p_thresh}, |FC| >= {args.fc_thresh}, "
             f"opposite-direction movement  |  AW-M-3 excluded  |  Chromatin: Keep+Review, SN: all proteins",
             ha="center", va="bottom", fontsize=6.5, color="#888888")

    # ---- source CSVs (membership per direction) ----
    def membership_table(sets):
        genes = sorted(set().union(*sets.values()))
        return pd.DataFrame([{"Gene": g, "AW": g in sets["AW"],
                              "Intox": g in sets["I"], "PA": g in sets["PA"]}
                             for g in genes])

    pthr = str(args.p_thresh).replace(".", "p")
    fcthr = str(args.fc_thresh).replace(".", "p")
    base = (f"{args.file_tag}_TrueTransloc_CrossConditionVenn_IntoChromAndSN"
            f"_p{pthr}_FC{fcthr}_{args.date}")
    out_pdf = os.path.join(args.outdir, base + ".pdf")
    out_png = os.path.join(args.outdir, base + ".png")
    membership_table(ch_sets).to_csv(
        os.path.join(args.outdir, base + "_IntoChromatin_membership.csv"), index=False)
    membership_table(sn_sets).to_csv(
        os.path.join(args.outdir, base + "_IntoSN_membership.csv"), index=False)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor="white")
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved:")
    print(" ", out_pdf)
    print(" ", out_png)


if __name__ == "__main__":
    main()
