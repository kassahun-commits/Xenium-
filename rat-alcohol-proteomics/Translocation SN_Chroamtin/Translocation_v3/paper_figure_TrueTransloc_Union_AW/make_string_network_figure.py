#!/usr/bin/env python3
"""
STRING protein-protein interaction network for a gene list (Rattus norvegicus),
drawn natively so text stays editable. One figure per direction (into chromatin
= up, into SN = down).

Style (cluster-coloured, after typical STRING / Cytoscape enrichment figures):
  * singletons (unconnected input proteins) are dropped;
  * every node is the SAME size dot;
  * nodes are COLOURED by functional module (greedy-modularity community);
  * only a subset of proteins are labelled (module hubs + small modules);
  * each major module is annotated with its top enriched term + FDR, pulled
    live from the STRING enrichment endpoint.

Interactions and enrichment are pulled live from the STRING REST API; edge
width encodes STRING combined score.

Standing rules (MEWS lab): editable vector text (pdf.fonttype=42); source tables
(edge list + node modules + module functions) saved alongside the figure; no
hardcoded paths.

Usage:
  python3 make_string_network_figure.py \
    --genelist current/..._IntoChromatin_..._genelist_....txt \
    --tag IntoChromatin --direction up --title "Into Chromatin (AW)" \
    --outdir current --species 10116 --required-score 400 --date 2026-06-08
"""
import os
import io
import time
import textwrap
import argparse
import datetime as _dt

import numpy as np
import pandas as pd
import requests
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

STRING_NET = "https://string-db.org/api/tsv/network"
STRING_ENRICH = "https://string-db.org/api/tsv/enrichment"

# qualitative palette for functional modules (largest module -> red, like STRING)
PALETTE = ["#E15759", "#59A14F", "#EDC948", "#4E79A7", "#B07AA1",
           "#F28E2B", "#76B7B2", "#FF9DA7", "#9C755F", "#86BCB6",
           "#D37295", "#A0CBE8"]
NEUTRAL = "#BAB0AC"   # extra/small modules beyond the palette


def fetch_string_enrichment(genes, species):
    """Top functional term for a gene set via STRING (GO Process preferred).
    Returns {term, fdr, category} or None on any failure / no hit."""
    if len(genes) < 3:
        return None
    data = {"identifiers": "\r".join(sorted(genes)), "species": species,
            "caller_identity": "mews_lab_translocation"}
    try:
        r = requests.post(STRING_ENRICH, data=data, timeout=60)
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text), sep="\t")
    except Exception:
        return None
    if df.empty or "category" not in df.columns:
        return None
    for cat in ("Process", "RCTM", "KEGG", "WikiPathways", "Function", "Component"):
        sub = df[df["category"] == cat].sort_values("fdr")
        if len(sub):
            row = sub.iloc[0]
            return {"term": str(row["description"]), "fdr": float(row["fdr"]),
                    "category": cat}
    return None


def fetch_string_edges(genes, species, required_score):
    """Return DataFrame of edges (gene_a, gene_b, score) among the input genes."""
    data = {
        "identifiers": "\r".join(genes),
        "species": species,
        "required_score": required_score,
        "caller_identity": "mews_lab_translocation",
    }
    r = requests.post(STRING_NET, data=data, timeout=60)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["gene_a", "gene_b", "score"])
    out = df[["preferredName_A", "preferredName_B", "score"]].rename(
        columns={"preferredName_A": "gene_a", "preferredName_B": "gene_b"})
    out = out.drop_duplicates(subset=["gene_a", "gene_b"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genelist", required=True)
    ap.add_argument("--tag", required=True, help="short label for filenames")
    ap.add_argument("--direction", choices=["up", "down"], required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--species", type=int, default=10116, help="rat=10116")
    ap.add_argument("--required-score", type=int, default=400,
                    help="STRING confidence 0-1000 (400=medium, 700=high)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.genelist) as f:
        genes = [ln.strip() for ln in f if ln.strip()]
    print(f"{args.tag}: {len(genes)} input genes")

    edges = fetch_string_edges(genes, args.species, args.required_score)
    print(f"  STRING edges (score>={args.required_score/1000:.2f}): {len(edges)}")

    G = nx.Graph()
    for _, e in edges.iterrows():
        G.add_edge(e["gene_a"], e["gene_b"], weight=float(e["score"]))
    n_connected = G.number_of_nodes()
    n_singletons = len(genes) - n_connected
    print(f"  connected nodes: {n_connected}, singletons: {n_singletons}")

    deg = dict(G.degree())

    # --- group proteins into functional modules (greedy modularity) ---
    communities = []
    if G.number_of_edges() > 0:
        communities = sorted(greedy_modularity_communities(G, weight="weight"),
                             key=len, reverse=True)
    node_group = {}
    for gi, comm in enumerate(communities):
        for nd in comm:
            node_group[nd] = gi
    group_color = {gi: (PALETTE[gi] if gi < len(PALETTE) else NEUTRAL)
                   for gi in range(len(communities))}
    print(f"  modules (size>=2): {len(communities)}")

    # --- layout: lay out each connected module on its own, then shelf-pack
    # the modules into a compact arrangement (singletons are dropped). ---
    pos = {}
    if G.number_of_edges() > 0:
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        boxes = []  # (nodes, local_pos, w, h)
        for comp in comps:
            H = G.subgraph(comp)
            m = H.number_of_nodes()
            nodes = list(H.nodes())
            if m == 2:
                w, h = 1.6, 0.7      # dyads -> tidy horizontal dumbbell
                local = {nodes[0]: (0.0, h / 2), nodes[1]: (w, h / 2)}
                boxes.append((nodes, local, w, h))
                continue
            p = nx.spring_layout(H, seed=args.seed,
                                 k=2.1 / np.sqrt(m), iterations=600, weight="weight")
            xs = np.array([p[u][0] for u in nodes])
            ys = np.array([p[u][1] for u in nodes])
            xr = (xs.max() - xs.min()) or 1.0
            yr = (ys.max() - ys.min()) or 1.0
            s = float(np.sqrt(m))
            local = {u: (((p[u][0] - xs.min()) / xr) * s,
                         ((p[u][1] - ys.min()) / yr) * s) for u in nodes}
            boxes.append((nodes, local, s, s))

        # if one module dominates (e.g. the SN giant component), give it its own
        # row so the small modules stack underneath instead of widening the page
        giant_frac = len(comps[0]) / max(n_connected, 1)
        if giant_frac > 0.5:
            target_w = boxes[0][2] * 1.02
        else:
            target_w = np.sqrt(sum(w * h for *_, w, h in boxes)) * 1.6
        pad = 0.5
        x0 = y_top = row_h = 0.0
        for nodes, local, w, h in boxes:
            if x0 > 0 and x0 + w > target_w:
                x0 = 0.0
                y_top -= (row_h + pad)
                row_h = 0.0
            for u in nodes:
                lx, ly = local[u]
                pos[u] = (x0 + lx, y_top - (h - ly))   # top-anchored, grows down
            x0 += w + pad
            row_h = max(row_h, h)

    # --- functional term per major module, placed in the L/R margins ---
    xs_all = [p[0] for p in pos.values()] or [0.0]
    ys_all = [p[1] for p in pos.values()] or [0.0]
    nx0, nx1 = min(xs_all), max(xs_all)
    ny0, ny1 = min(ys_all), max(ys_all)
    gcx = 0.5 * (nx0 + nx1)
    IPU = 1.0  # inches per data unit

    raw = []
    for gi, comm in enumerate(communities):
        if len(comm) < 4:
            continue
        enr = fetch_string_enrichment(list(comm), args.species)
        time.sleep(0.4)            # be gentle to the STRING API
        if not enr:
            continue
        cx = float(np.mean([pos[n][0] for n in comm]))
        cy = float(np.mean([pos[n][1] for n in comm]))
        raw.append({"group": gi, "size": len(comm),
                    "term": enr["term"], "term_disp": "\n".join(textwrap.wrap(
                        enr["term"], 24)) or enr["term"],
                    "fdr": enr["fdr"], "category": enr["category"],
                    "cx": cx, "cy": cy, "color": group_color[gi]})
        print(f"  module {gi} (n={len(comm)}): {enr['term']} (FDR={enr['fdr']:.1e})")

    # assign each module label to the nearer margin, spread out vertically
    span_w = (nx1 - nx0) or 1.0
    gap = 0.13 * span_w
    left = sorted([c for c in raw if c["cx"] < gcx], key=lambda c: -c["cy"])
    right = sorted([c for c in raw if c["cx"] >= gcx], key=lambda c: -c["cy"])
    cluster_funcs = []
    for items, lx, ha in ((left, nx0 - gap, "right"), (right, nx1 + gap, "left")):
        if not items:
            continue
        slots = np.linspace(ny1, ny0, len(items)) if len(items) > 1 \
            else [0.5 * (ny0 + ny1)]
        for c, ly in zip(items, slots):
            c.update(lx=lx, ly=float(ly), ha=ha)
            cluster_funcs.append(c)

    # --- size the figure to include nodes + (wrapped) margin labels ---
    def _tw(disp):
        longest = max([len(s) for s in disp.split("\n")] + [len("FDR=0.0e-00")])
        return longest * 11.5 * 0.62 / 72.0 / IPU   # ~text width in data units
    fx, fy = [nx0, nx1], [ny0, ny1]
    for c in cluster_funcs:
        w = _tw(c["term_disp"])
        fx += [c["lx"], c["lx"] + (w if c["ha"] == "left" else -w)]
        fy += [c["ly"] + 0.5, c["ly"] - 0.5]
    pad_d = 0.04 * max((max(fx) - min(fx)), (max(fy) - min(fy)))
    xmin, xmax = min(fx) - pad_d, max(fx) + pad_d
    ymin, ymax = min(fy) - pad_d, max(fy) + pad_d
    base_w = float(np.clip((xmax - xmin) * IPU, 9.0, 46.0))
    fig_h = float(np.clip((ymax - ymin) * IPU, 6.0, 36.0)) + 1.0
    fig = plt.figure(figsize=(base_w, fig_h), facecolor="white")
    ax = fig.add_axes([0.02, 0.01, 0.96, 1.0 - 1.0 / fig_h])
    ax.axis("off")

    if G.number_of_edges() > 0:
        node_s = float(np.clip(260 - 0.5 * n_connected, 110, 240))
        # --- edges (grey, thicker = higher STRING confidence) ---
        for (u, v) in G.edges():
            s = G[u][v]["weight"]
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color="#9AA0A6", lw=0.4 + 1.8 * (s - 0.4) / 0.6,
                    alpha=0.35 + 0.4 * (s - 0.4) / 0.6, zorder=1)
        # --- nodes: uniform size, coloured by functional module ---
        ncols = [group_color.get(node_group.get(x), NEUTRAL) for x in G.nodes()]
        ax.scatter([pos[x][0] for x in G.nodes()], [pos[x][1] for x in G.nodes()],
                   s=node_s, c=ncols, edgecolors="white", linewidths=1.0, zorder=2)
        # --- protein labels: module hubs + members of tiny modules only ---
        comp_size = {}
        for comp in nx.connected_components(G):
            for nd in comp:
                comp_size[nd] = len(comp)
        to_label = set()
        for comm in communities:
            members = sorted(comm, key=lambda d: deg[d], reverse=True)
            k = max(2, int(round(0.35 * len(members))))
            to_label.update(members[:k])
        to_label.update([nd for nd in G.nodes() if comp_size[nd] <= 3])
        hubs = {sorted(comm, key=lambda d: deg[d])[-1] for comm in communities if comm}
        for x in to_label:
            ax.text(pos[x][0], pos[x][1] + 0.07, x,
                    fontsize=9.0 if x in hubs else 7.2,
                    fontweight="bold" if x in hubs else "normal",
                    ha="center", va="bottom", zorder=4,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
        # --- functional module labels in the margins, leader line to module ---
        for c in cluster_funcs:
            ax.annotate(f"{c['term_disp']}\nFDR={c['fdr']:.1e}",
                        xy=(c["cx"], c["cy"]), xytext=(c["lx"], c["ly"]),
                        ha=c["ha"], va="center", fontsize=11.5, color=c["color"],
                        fontweight="bold", zorder=5,
                        path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
                        arrowprops=dict(arrowstyle="-", color=c["color"], lw=0.9,
                                        alpha=0.55, shrinkA=4, shrinkB=4))
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")

    # title + stats
    fig.text(0.5, 0.985, f"{args.title} — STRING interaction network",
             ha="center", va="top", fontsize=16, fontweight="bold", color="#1A3A5C")
    fig.text(0.5, 0.962,
             f"{n_connected} connected proteins in {len(communities)} modules  |  "
             f"{G.number_of_edges()} edges  |  STRING score ≥ {args.required_score/1000:.2f}"
             f"  |  {n_singletons} singletons not shown",
             ha="center", va="top", fontsize=9, color="#666666")

    base = (f"{args.file_tag}_AW_{args.tag}_STRINGnetwork"
            f"_score{args.required_score}_{args.date}")
    out_pdf = os.path.join(args.outdir, base + ".pdf")
    out_png = os.path.join(args.outdir, base + ".png")
    edges.to_csv(os.path.join(args.outdir, base + "_edges.csv"), index=False)
    pd.DataFrame({"gene": list(deg.keys()),
                  "degree": [deg[g] for g in deg],
                  "module": [node_group.get(g, -1) for g in deg]}) \
        .sort_values(["module", "degree"], ascending=[True, False]) \
        .to_csv(os.path.join(args.outdir, base + "_node_modules.csv"), index=False)
    pd.DataFrame(cluster_funcs,
                 columns=["group", "size", "category", "term", "fdr"]) \
        .to_csv(os.path.join(args.outdir, base + "_module_functions.csv"), index=False)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor="white")
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved:")
    print(" ", out_pdf)
    print(" ", out_png)


if __name__ == "__main__":
    main()
