#!/usr/bin/env python3
"""
STRING protein-protein interaction network, "flat" Cytoscape/STRING-panel style
(Rattus norvegicus), drawn natively so all text stays editable.

Design goal: match a compact published STRING panel (uniform-colour nodes, soft
darker node ring, grey edges whose width + darkness encode the STRING combined
score, every node labelled next to its dot, no module colouring, white
background). One figure per direction:
  * up   = into chromatin  -> all nodes RED
  * down = into SN         -> all nodes BLUE

Layout: the largest connected component is laid out with a force-directed
(spring) layout; the remaining small components (dyads/triads) are ring-packed
tightly around it so the whole thing reads as one compact cluster like the
reference panel. Singletons (unconnected input proteins) are dropped.

Edges/scores can be re-fetched live from the STRING REST API, or (default when
--edges is given) read from the previously saved edge list so the network is
byte-for-byte identical to the earlier figure - only the styling changes.

Standing rules (MEWS lab): editable vector text (pdf.fonttype=42, ps.fonttype=42,
Arial); source tables (edge list + node list) saved alongside the figure; no
hardcoded paths; current/ + previous_versions/ structure handled by the caller.

Usage:
  python3 make_string_network_flat.py \
    --genelist current/..._IntoChromatin_..._genelist_....txt \
    --edges    current/..._IntoChromatin_..._STRINGnetwork_..._edges.csv \
    --tag IntoChromatin --direction up --title "Into Chromatin (AW)" \
    --outdir current --species 10116 --required-score 400 --date 2026-06-08
"""
import os
import argparse
import datetime as _dt

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"

# ---- style block (single source of truth for the figure look) ----------------
NODE_STYLE = {
    "up":   {"face": "#E15759", "edge": "#9C2B2A"},   # into chromatin -> red
    "down": {"face": "#4E79A7", "edge": "#2C4B6E"},   # into SN        -> blue
}
EDGE_LIGHT = np.array([0.82, 0.82, 0.83])   # low-confidence edge colour
EDGE_DARK = np.array([0.27, 0.27, 0.34])    # high-confidence edge colour
EDGE_LW = (0.5, 3.2)                          # linewidth at score 0.4 .. 1.0
NODE_RING_LW = 1.3
LABEL_COLOR = "#1A1A1A"
LABEL_HALO_LW = 2.2
SCORE_MIN = 0.4                               # required_score / 1000


def load_edges(args):
    """Edge DataFrame (gene_a, gene_b, score in 0..1) from CSV or live STRING."""
    if args.edges:
        df = pd.read_csv(args.edges)
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        df = df.dropna(subset=["score"])[["gene_a", "gene_b", "score"]]
    else:
        # live fall-back: reuse the fetcher from the original module
        from make_string_network_figure import fetch_string_edges
        with open(args.genelist) as f:
            genes = [ln.strip() for ln in f if ln.strip()]
        df = fetch_string_edges(genes, args.species, args.required_score)
        df["score"] = pd.to_numeric(df["score"], errors="coerce") / 1000.0
        df = df.dropna(subset=["score"])
    if args.min_score:
        df = df[df["score"] >= args.min_score]
    return df.reset_index(drop=True)


def _layout_component(H, seed):
    """Local layout for one component, min-corner at (0,0); returns (pos, w, h)."""
    nodes = list(H.nodes())
    m = len(nodes)
    if m == 1:
        return {nodes[0]: np.array([0.0, 0.0])}, 0.0, 0.0
    if m == 2:
        d = 1.0
        p = {nodes[0]: np.array([0.0, 0.0]), nodes[1]: np.array([d, 0.0])}
        return p, d, 0.0
    # more internal spread for big components so labels have room
    raw = nx.spring_layout(H, seed=seed, k=2.6 / np.sqrt(m),
                           iterations=700, weight="weight")
    pts = np.array([raw[u] for u in nodes])
    span = np.ptp(pts, axis=0)
    span[span == 0] = 1.0
    pts = pts / span.max() * (1.25 * np.sqrt(m))   # uniform inter-node spacing
    pts -= pts.min(axis=0)                          # min corner -> origin
    w, h = float(pts[:, 0].max()), float(pts[:, 1].max())
    return {u: pts[i] for i, u in enumerate(nodes)}, w, h


def _skyline_pack(items, W, gap):
    """Bottom-left skyline packing of items [(lp, w, h)] into a strip of width
    W; components drop into the lowest gap that fits, so the result is a tight
    rectangle with the empty space filled. Returns (pos, used_w, used_h)."""
    sky = [[0.0, 0.0]]   # breakpoints [x, height]; each spans to the next x

    def height_at(x):
        h = sky[0][1]
        for sx, sh in sky:
            if sx <= x + 1e-9:
                h = sh
            else:
                break
        return h

    def max_height(xl, xr):
        h = height_at(xl)
        for sx, sh in sky:
            if xl < sx < xr:
                h = max(h, sh)
        return h

    def raise_region(xl, xr, h):
        h_after = height_at(xr)
        keep = [n for n in sky if not (xl - 1e-9 < n[0] < xr + 1e-9)]
        keep.append([xl, h])
        keep.append([xr, h_after])
        keep.sort()
        out = []
        for n in keep:
            if out and abs(out[-1][0] - n[0]) < 1e-9:
                out[-1] = n
            else:
                out.append(n)
        sky[:] = out

    pos = {}
    used_w = used_h = 0.0
    for lp, w, h in items:
        wpad = w + gap
        cands = sorted({n[0] for n in sky if n[0] + wpad <= W + 1e-9}) or [0.0]
        best_x = best_y = None
        for cx in cands:
            y = max_height(cx, cx + wpad)
            if best_y is None or y < best_y - 1e-9 or \
                    (abs(y - best_y) < 1e-9 and cx < best_x):
                best_x, best_y = cx, y
        for u, p in lp.items():
            pos[u] = np.array([best_x + p[0], best_y + p[1]])
        raise_region(best_x, best_x + wpad, best_y + h + gap)
        used_w = max(used_w, best_x + w)
        used_h = max(used_h, best_y + h)
    return pos, used_w, used_h


def compact_layout(G, seed):
    """Lay out every connected component, then skyline-pack them all into one
    tight rectangle (near-square aspect), so the whole network sits on a single
    'page' with minimal whitespace regardless of how fragmented it is."""
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    boxes = []
    for c in comps:
        lp, w, h = _layout_component(G.subgraph(c), seed)
        boxes.append((lp, max(w, 0.6), max(h, 0.6)))
    boxes.sort(key=lambda b: b[2], reverse=True)        # tallest first

    gap = 0.7
    area = sum((w + gap) * (h + gap) for _, w, h in boxes)
    min_w = max(w for _, w, _ in boxes)
    best = None
    for f in (0.9, 1.1, 1.3, 1.6, 1.9, 2.3):            # try several widths
        W = max(min_w, float(np.sqrt(area * f)))
        pos, uw, uh = _skyline_pack(boxes, W, gap)
        score = (uw * uh) * (1.0 + 0.12 * abs(uw / max(uh, 1e-6) - 1.35))
        if best is None or score < best[0]:
            best = (score, pos)
    return best[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genelist", required=True)
    ap.add_argument("--edges", default=None,
                    help="saved edge list (gene_a,gene_b,score 0..1); else fetch live")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--direction", choices=["up", "down"], required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--species", type=int, default=10116)
    ap.add_argument("--required-score", type=int, default=400)
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="extra edge filter on the loaded list (e.g. 0.7=high conf)")
    ap.add_argument("--largest-only", action="store_true",
                    help="keep only the largest connected component (clean blob)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--date", default=_dt.date.today().isoformat())
    ap.add_argument("--file-tag", default="RatAlcoholProteome")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.genelist) as f:
        genes = [ln.strip() for ln in f if ln.strip()]
    edges = load_edges(args)

    G = nx.Graph()
    for _, e in edges.iterrows():
        G.add_edge(e["gene_a"], e["gene_b"], weight=float(e["score"]))
    if args.largest_only and G.number_of_nodes():
        giant = max(nx.connected_components(G), key=len)
        G = G.subgraph(giant).copy()
        edges = edges[edges["gene_a"].isin(giant) & edges["gene_b"].isin(giant)]
    n_conn = G.number_of_nodes()
    n_single = len(set(genes)) - n_conn
    print(f"{args.tag}: {len(genes)} input | {n_conn} connected | "
          f"{G.number_of_edges()} edges | {n_single} singletons dropped")

    pos = compact_layout(G, args.seed)
    deg = dict(G.degree())

    # --- figure size from layout extent (compact, aspect-equal) ---------------
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    pad = 1.2
    xmin, xmax = xs.min() - pad, xs.max() + pad
    ymin, ymax = ys.min() - pad, ys.max() + pad
    span = max(xmax - xmin, ymax - ymin)
    in_per_unit = float(np.clip(11.0 / np.sqrt(max(n_conn, 1)), 0.42, 1.0))
    fig_w = float(np.clip((xmax - xmin) * in_per_unit, 8.0, 26.0))
    fig_h = float(np.clip((ymax - ymin) * in_per_unit, 6.0, 26.0)) + 0.7
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    ax = fig.add_axes([0.01, 0.01, 0.98, 1.0 - 1.0 / fig_h])
    ax.axis("off")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")

    # --- edges: grey, width + darkness scale with STRING combined score -------
    for u, v in G.edges():
        s = G[u][v]["weight"]
        t = float(np.clip((s - SCORE_MIN) / (1.0 - SCORE_MIN), 0.0, 1.0))
        col = EDGE_LIGHT + t * (EDGE_DARK - EDGE_LIGHT)
        lw = EDGE_LW[0] + t * (EDGE_LW[1] - EDGE_LW[0])
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=col, lw=lw, alpha=0.45 + 0.45 * t, zorder=1,
                solid_capstyle="round")

    # --- nodes: uniform colour per direction, soft darker ring ----------------
    st = NODE_STYLE[args.direction]
    node_s = float(np.clip(1900.0 / np.sqrt(max(n_conn, 1)), 70.0, 470.0))
    ax.scatter([pos[x][0] for x in G.nodes()], [pos[x][1] for x in G.nodes()],
               s=node_s, c=st["face"], edgecolors=st["edge"],
               linewidths=NODE_RING_LW, zorder=2)

    # --- labels: every node, glued just above its OWN circle, white halo ------
    # convert the node radius (points) into data units so the label sits exactly
    # on the rim of its dot at any figure scale (fixes "floating" labels).
    ax_w_in, ax_h_in = fig_w * 0.98, fig_h - 1.0
    data_per_pt = max((xmax - xmin) / (ax_w_in * 72.0),
                      (ymax - ymin) / (ax_h_in * 72.0))
    node_r_data = np.sqrt(node_s / np.pi) * data_per_pt
    off = node_r_data * 1.15
    lab_fs = float(np.clip(150.0 / np.sqrt(max(n_conn, 1)), 5.0, 10.0))
    for x in G.nodes():
        ax.text(pos[x][0], pos[x][1] + off, x, fontsize=lab_fs,
                ha="center", va="bottom", color=LABEL_COLOR, zorder=4,
                path_effects=[pe.withStroke(linewidth=LABEL_HALO_LW,
                                            foreground="white")])

    # title + one-line stats
    fig.text(0.5, 0.985, f"{args.title} — STRING interaction network",
             ha="center", va="top", fontsize=15, fontweight="bold",
             color="#1A3A5C")
    fig.text(0.5, 0.96,
             f"{n_conn} connected proteins  |  {G.number_of_edges()} edges  |  "
             f"STRING score ≥ {args.required_score/1000:.2f}  |  "
             f"{n_single} singletons not shown",
             ha="center", va="top", fontsize=9, color="#666666")

    base = (f"{args.file_tag}_AW_{args.tag}_STRINGnetwork"
            f"_score{args.required_score}_{args.date}")
    out_pdf = os.path.join(args.outdir, base + ".pdf")
    out_png = os.path.join(args.outdir, base + ".png")
    edges.to_csv(os.path.join(args.outdir, base + "_edges.csv"), index=False)
    pd.DataFrame({"gene": list(G.nodes()),
                  "degree": [deg[g] for g in G.nodes()]}) \
        .sort_values("degree", ascending=False) \
        .to_csv(os.path.join(args.outdir, base + "_nodes.csv"), index=False)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, facecolor="white")
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved:\n ", out_pdf, "\n ", out_png)


if __name__ == "__main__":
    main()
