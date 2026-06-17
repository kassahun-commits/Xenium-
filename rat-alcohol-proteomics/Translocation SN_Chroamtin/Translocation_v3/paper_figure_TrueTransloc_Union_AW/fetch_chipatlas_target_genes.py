#!/usr/bin/env python3
"""
ChIP-Atlas Target Genes for AW true translocators.

For each translocating protein that has public ChIP-seq experiments in
ChIP-Atlas, pull the precomputed "Target Genes" (genes with binding peaks
within +/- <distance> kb of their TSS) and rank them by the averaged binding
score across all experiments.

Rationale: proteins LEAVING chromatin (into Soluble Nuclear) during acute
withdrawal are enriched for chromatin/transcriptional regulators (Creb1, Mta1,
Smarcc2, Carm1, Cbx1, Top2a, ...). Their canonical target genes give a
mechanistic hypothesis for which downstream genes are affected when they
redistribute.

Caveats (reported, not hidden): ChIP-Atlas data is overwhelmingly human/mouse
cell lines (not rat brain), so these are CANONICAL targets, not AW-specific.
Rat (rn6) coverage is sparse, so we query human (hg38) first, then mouse (mm10),
then rat (rn6).

Data source (static files, no key needed):
  discovery : https://chip-atlas.dbcls.jp/data/metadata/analysisList.tab
  targets   : https://chip-atlas.dbcls.jp/data/<genome>/target/<ANTIGEN>.<dist>.tsv
              (col0 = Target_genes, col1 = '<ANTIGEN>|Average', then one column
               per experiment; values are averaged MACS2 binding scores)

Standing rules (MEWS lab): source CSVs alongside outputs, no hardcoded paths,
current/ + previous_versions/ structure.

Usage:
  python3 fetch_chipatlas_target_genes.py \
    --genelist up=current/..._IntoChromatin_..._genelist_....txt \
    --genelist down=current/..._IntoSN_..._genelist_....txt \
    --outdir ChIPAtlas_TargetGenes/current --distance 5 --top-n 30 \
    --date 2026-06-08
"""
import os
import sys
import argparse
import datetime as _dt
import urllib.request

import numpy as np
import pandas as pd

CHIP_BASE = "https://chip-atlas.dbcls.jp/data"
ANALYSIS_LIST = f"{CHIP_BASE}/metadata/analysisList.tab"
GENOME_PRIORITY = ["hg38", "mm10", "rn6", "hg19", "mm9"]


def load_analysis_sets(cache):
    path = os.path.join(cache, "analysisList.tab")
    if not os.path.exists(path):
        print("Downloading analysisList.tab ...")
        urllib.request.urlretrieve(ANALYSIS_LIST, path)
    sets = {}
    with open(path) as f:
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 4:
                sets.setdefault(c[3], set()).add(c[0])
    return sets


def match_genome(gene, sets):
    """Return (genome, antigen) for the highest-priority genome that has data."""
    upper = gene.upper()
    title = gene[:1].upper() + gene[1:].lower()
    cand = {"hg38": upper, "hg19": upper,
            "mm10": title, "mm9": title, "rn6": title}
    for g in GENOME_PRIORITY:
        ag = cand.get(g)
        if ag and ag in sets.get(g, set()):
            return g, ag
    return None, None


def fetch_target_file(genome, antigen, dist, cache):
    fn = f"{antigen}.{dist}.tsv"
    path = os.path.join(cache, f"{genome}_{fn}")
    if not os.path.exists(path):
        url = f"{CHIP_BASE}/{genome}/target/{fn}"
        urllib.request.urlretrieve(url, path)
    return path


def summarize_targets(path):
    """Return (ranked_df[gene,avg_score], n_experiments)."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    n_exp = max(0, len(header) - 2)  # minus Target_genes + Average
    df = pd.read_csv(path, sep="\t", usecols=[0, 1], dtype={0: str})
    df.columns = ["target_gene", "avg_score"]
    df["avg_score"] = pd.to_numeric(df["avg_score"], errors="coerce")
    df = df.dropna(subset=["avg_score"])
    df = df[df["avg_score"] > 0].sort_values("avg_score", ascending=False)
    df = df.reset_index(drop=True)
    return df, n_exp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genelist", action="append", required=True,
                    help="label=path/to/genelist.txt  (repeatable, e.g. up=..., down=...)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--distance", type=int, default=5, choices=[1, 5, 10],
                    help="kb from TSS for peak-to-gene assignment")
    ap.add_argument("--top-n", type=int, default=30,
                    help="number of top target genes to name in the summary")
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cache = os.path.join(args.outdir, "_raw_targetgenes")
    os.makedirs(cache, exist_ok=True)

    sets = load_analysis_sets(cache)

    for spec in args.genelist:
        if "=" not in spec:
            sys.exit(f"--genelist must be label=path, got: {spec}")
        label, path = spec.split("=", 1)
        genes = [l.strip() for l in open(path) if l.strip()]
        print(f"\n===== {label}: {len(genes)} input proteins =====")

        rows = []
        union_targets = {}  # target_gene -> set(regulators)
        for g in genes:
            genome, antigen = match_genome(g, sets)
            if genome is None:
                continue
            try:
                tpath = fetch_target_file(genome, antigen, args.distance, cache)
                tdf, n_exp = summarize_targets(tpath)
            except Exception as e:
                print(f"  ! {g} ({genome}:{antigen}) fetch failed: {e}")
                continue
            top = tdf.head(args.top_n)["target_gene"].tolist()
            rows.append({
                "protein": g, "genome": genome, "antigen": antigen,
                "n_experiments": n_exp, "n_target_genes": len(tdf),
                f"top{args.top_n}_target_genes": "; ".join(top),
            })
            # per-protein full ranked CSV
            per = os.path.join(
                args.outdir,
                f"{args.file_tag}_ChIPAtlas_{label}_{g}_{genome}_d{args.distance}kb_targetgenes.csv")
            tdf.to_csv(per, index=False)
            for tg in top:
                union_targets.setdefault(tg.upper(), set()).add(g)
            print(f"  {g:9s} {genome}:{antigen:9s}  exps={n_exp:3d}  "
                  f"targets(avg>0)={len(tdf):5d}  top: {', '.join(top[:6])}")

        if not rows:
            print(f"  (no proteins in '{label}' had ChIP-Atlas target-gene data)")
            continue

        # per-direction summary
        summ = pd.DataFrame(rows).sort_values("n_experiments", ascending=False)
        base = os.path.join(
            args.outdir,
            f"{args.file_tag}_ChIPAtlas_{label}_TargetGenesSummary_d{args.distance}kb_{args.date}")
        summ.to_csv(base + ".csv", index=False)

        # union target genes with how many regulators hit them (recurrent targets)
        ut = pd.DataFrame(
            [{"target_gene": k, "n_regulators": len(v),
              "regulators": "; ".join(sorted(v))} for k, v in union_targets.items()]
        ).sort_values(["n_regulators", "target_gene"], ascending=[False, True])
        ut.to_csv(base + "_unionTargets.csv", index=False)

        print(f"  -> {len(rows)} proteins with data; "
              f"{len(ut)} unique top-{args.top_n} target genes "
              f"(max regulators on one gene = {int(ut['n_regulators'].max())})")
        print(f"  wrote: {base}.csv")
        print(f"         {base}_unionTargets.csv")


if __name__ == "__main__":
    main()
