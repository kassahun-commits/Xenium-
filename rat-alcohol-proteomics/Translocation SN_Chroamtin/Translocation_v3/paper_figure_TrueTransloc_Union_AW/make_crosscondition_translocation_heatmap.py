#!/usr/bin/env python3
"""
Cross-condition translocation heatmap (Intoxication -> Acute Withdrawal ->
Protracted Abstinence) of TRUE Chromatin <-> Soluble-Nuclear translocators.

Motivation
----------
The cross-condition Venn collapses two things we care about: the *magnitude*
of translocation and *direction switching*. I / AW / PA form a disease-course
trajectory, so this script visualizes, for every protein that is a true
translocator in >=1 condition, the signed translocation effect size across the
three conditions:

    cell value = Interaction = mean(delta_chromatin) - mean(delta_SN)
                 (positive -> net shift INTO chromatin  [red],
                  negative -> net shift INTO soluble nuclear [blue])

Significance (BH FDR p_adj < threshold) is overlaid as a dot; rows are
hierarchically clustered; a left color strip classifies each protein as:
    Transient (true in exactly one condition)
    Shared    (true, same direction, in >=2 conditions)
    Switch    (true INTO chromatin in one condition and INTO SN in another)

The statistics engine (build_cond_stats) is imported verbatim from
make_crosscondition_translocation_venns.py so thresholds, the AW-M-3 exclusion,
Chromatin Keep+Review filtering and SN-all handling stay identical to the
Venns / Panels A-B.

Standing rules (MEWS lab): editable vector text (pdf.fonttype=42), source CSVs
alongside the figure, no hardcoded paths.

Usage:
  python3 make_crosscondition_translocation_heatmap.py \
    --xlsx "/.../EDIT  Excluding AWM3 ... copy.xlsx" \
    --outdir current --p-thresh 0.10 --fc-thresh 0.5 --date 2026-06-08
"""
import os
import argparse
import datetime as _dt

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list

# identical statistics engine as the Venns / Panels (single source of truth)
from make_crosscondition_translocation_venns import build_cond_stats

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

CONDS = [("I", "Intox."), ("AW", "Acute\nWithdrawal"), ("PA", "Protracted\nAbstinence")]
CLASS_ORDER = ["Switch", "Shared", "Transient"]
CLASS_COLORS = {"Transient": "#D9D9D9", "Shared": "#59A14F", "Switch": "#F28E2B"}
CMAP = "RdBu_r"   # red = into chromatin (+), blue = into SN (-)


def classify_row(dirs):
    """dirs = list of TrueDir strings across conditions (order I, AW, PA)."""
    present = [d for d in dirs if d in ("Into_Chromatin", "Into_SN")]
    if len(set(present)) == 2:
        return "Switch"
    if len(present) >= 2:
        return "Shared"
    return "Transient"


def build_matrix(ch_df, sn_df, common, p_thresh, fc_thresh):
    inter, padj, tdir = {}, {}, {}
    for c, _ in CONDS:
        d = build_cond_stats(ch_df, sn_df, common, c, p_thresh, fc_thresh).set_index("Gene")
        inter[c], padj[c], tdir[c] = d["Interaction"], d["p_adj"], d["TrueDir"]
        nch = int((d["TrueDir"] == "Into_Chromatin").sum())
        nsn = int((d["TrueDir"] == "Into_SN").sum())
        print(f"  {c}: Into Chromatin={nch}, Into SN={nsn}")

    is_true = pd.DataFrame({c: tdir[c] != "NS" for c, _ in CONDS})
    union = is_true.index[is_true.any(axis=1)].tolist()
    print(f"Union of true translocators across I/AW/PA: {len(union)}")

    M = pd.DataFrame({c: inter[c].reindex(union) for c, _ in CONDS})
    P = pd.DataFrame({c: padj[c].reindex(union) for c, _ in CONDS})
    Dr = pd.DataFrame({c: tdir[c].reindex(union) for c, _ in CONDS})
    M.columns = P.columns = Dr.columns = [c for c, _ in CONDS]

    cls = Dr.apply(lambda r: classify_row(list(r.values)), axis=1)
    n_true = (Dr != "NS").sum(axis=1)
    # AW direction (for stable secondary sort / reporting)
    aw_dir = Dr["AW"].where(Dr["AW"] != "NS", "")
    meta = pd.DataFrame({"Class": cls, "n_conditions_true": n_true, "AW_dir": aw_dir})
    return M, P, Dr, meta


def cluster_order(M):
    X = M.fillna(0.0).values
    if X.shape[0] < 3:
        return np.arange(X.shape[0])
    Z = linkage(X, method="ward")
    return leaves_list(Z)


def draw_heatmap(M, P, meta, p_thresh, title, outbase, label_genes):
    n, ncol = M.shape
    genes = M.index.tolist()
    vals = np.ma.masked_invalid(M.values)
    vlim = np.nanpercentile(np.abs(M.values), 98)
    vlim = float(vlim) if np.isfinite(vlim) and vlim > 0 else 1.0
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0.0, vmax=vlim)
    cmap = plt.get_cmap(CMAP).copy()
    cmap.set_bad("#EDEDED")

    row_h = 0.15 if label_genes else 0.045
    fig_h = max(3.2, row_h * n + 1.7)
    fig = plt.figure(figsize=(5.6, fig_h), facecolor="white")

    left = 0.40 if label_genes else 0.06
    strip_w = 0.05
    hm_w = 0.34
    ax_strip = fig.add_axes([left, 0.10, strip_w, 0.80])
    ax = fig.add_axes([left + strip_w + 0.012, 0.10, hm_w, 0.80])
    ax_cb = fig.add_axes([left + strip_w + 0.012 + hm_w + 0.03, 0.10, 0.025, 0.80])

    # ---- class strip ----
    code = {cl: i for i, cl in enumerate(CLASS_ORDER)}
    strip_codes = np.array([[code[c]] for c in meta["Class"].values])
    strip_cmap = ListedColormap([CLASS_COLORS[c] for c in CLASS_ORDER])
    ax_strip.imshow(strip_codes, aspect="auto", cmap=strip_cmap,
                    vmin=-0.5, vmax=len(CLASS_ORDER) - 0.5)
    ax_strip.set_xticks([])
    ax_strip.set_yticks([])
    ax_strip.set_title("class", fontsize=7, pad=4)

    # ---- heatmap ----
    im = ax.imshow(vals, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels([lbl for _, lbl in CONDS], fontsize=9)
    ax.set_yticks([])
    # significance dots
    ys, xs = np.where((P.values < p_thresh))
    ax.scatter(xs, ys, s=6 if not label_genes else 11, c="#111111",
               edgecolors="white", linewidths=0.3, zorder=4)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if label_genes:
        ax.set_yticks(range(n))
        ax.set_yticklabels(genes, fontsize=6.2)
        ax.tick_params(axis="y", length=0)

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("Translocation score\n(+ into chromatin / − into SN)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    fig.suptitle(title, fontsize=11, fontweight="bold", color="#1A3A5C", y=0.985)

    # legend for class strip + significance
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", markersize=8,
                          markerfacecolor=CLASS_COLORS[c], markeredgecolor="none",
                          label=c) for c in CLASS_ORDER]
    handles.append(plt.Line2D([0], [0], marker="o", linestyle="", markersize=6,
                              markerfacecolor="#111111", markeredgecolor="white",
                              label=f"p_adj < {p_thresh}"))
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, 0.045))

    with PdfPages(outbase + ".pdf") as pdf:
        pdf.savefig(fig, facecolor="white", bbox_inches="tight")
    fig.savefig(outbase + ".png", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def draw_class_bar(meta, outbase):
    # counts per class, split by AW direction where available
    fig, ax = plt.subplots(figsize=(4.6, 3.2), facecolor="white")
    counts = meta["Class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    bars = ax.bar(range(len(CLASS_ORDER)), counts.values,
                  color=[CLASS_COLORS[c] for c in CLASS_ORDER], zorder=3)
    for b, v in zip(bars, counts.values):
        ax.text(b.get_x() + b.get_width() / 2, v, str(v), ha="center", va="bottom",
                fontsize=10)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels(CLASS_ORDER, fontsize=10)
    ax.set_ylabel("True translocators (proteins)", fontsize=10)
    ax.set_title("Cross-condition translocation persistence", fontsize=11,
                 fontweight="bold", color="#1A3A5C")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, counts.max() * 1.15 if counts.max() > 0 else 1)
    fig.tight_layout()
    with PdfPages(outbase + ".pdf") as pdf:
        pdf.savefig(fig, facecolor="white")
    fig.savefig(outbase + ".png", dpi=300, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--p-thresh", type=float, default=0.10)
    ap.add_argument("--fc-thresh", type=float, default=0.5)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    ap.add_argument("--label-max", type=int, default=90,
                    help="label every gene on the full heatmap if union <= this")
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

    M, P, Dr, meta = build_matrix(ch_df, sn_df, common, args.p_thresh, args.fc_thresh)

    # cluster full union
    order = cluster_order(M)
    M, P, Dr, meta = M.iloc[order], P.iloc[order], Dr.iloc[order], meta.iloc[order]

    pthr = str(args.p_thresh).replace(".", "p")
    fcthr = str(args.fc_thresh).replace(".", "p")
    base = os.path.join(
        args.outdir,
        f"{args.file_tag}_TrueTransloc_CrossConditionHeatmap_IAWvsPA"
        f"_p{pthr}_FC{fcthr}_{args.date}")

    # ---- full landscape heatmap ----
    n = M.shape[0]
    draw_heatmap(M, P, meta, args.p_thresh,
                 f"Cross-condition translocation (n={n} true translocators)",
                 base + "_full", label_genes=(n <= args.label_max))

    # ---- labeled cross-condition core (Shared + Switch) ----
    core_mask = meta["Class"].isin(["Shared", "Switch"])
    if core_mask.sum() >= 2:
        Mc, Pc, mc = M[core_mask], P[core_mask], meta[core_mask]
        oc = cluster_order(Mc)
        Mc, Pc, mc = Mc.iloc[oc], Pc.iloc[oc], mc.iloc[oc]
        draw_heatmap(Mc, Pc, mc, args.p_thresh,
                     f"Cross-condition core: Shared + Switch (n={Mc.shape[0]})",
                     base + "_core", label_genes=True)
        print(f"Core (Shared+Switch): {int(core_mask.sum())} proteins")
    else:
        print("Core (Shared+Switch): <2 proteins, core heatmap skipped")

    # ---- persistence bar ----
    draw_class_bar(meta, base + "_classcounts")

    # ---- source CSVs ----
    src = pd.DataFrame({"Gene": M.index})
    for c, _ in CONDS:
        src[f"{c}_Interaction"] = M[c].values
        src[f"{c}_p_adj"] = P[c].values
        src[f"{c}_TrueDir"] = Dr[c].values
    src["n_conditions_true"] = meta["n_conditions_true"].values
    src["Class"] = meta["Class"].values
    src["AW_dir"] = meta["AW_dir"].values
    src.to_csv(base + "_source.csv", index=False)
    meta["Class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int) \
        .rename_axis("Class").rename("count").to_csv(base + "_classcounts_source.csv")

    print("\nSaved:")
    for suff in ("_full.pdf", "_core.pdf", "_classcounts.pdf", "_source.csv"):
        print("  ", base + suff)


if __name__ == "__main__":
    main()
