#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled — cell-type UMAP + marker / MCT-axis feature grid.

Grant figure for the neurons-vs-astrocytes story. Two pages per version:

  PAGE A  cell-type UMAP — one embedding coloured by cell type
          (Excitatory / Inhibitory / Astrocyte / Oligodendrocyte / Microglia /
          Other). Shows that neurons and astrocytes occupy transcriptionally
          distinct territory.
  PAGE B  the SAME embedding shown as a 3x4 grid of feature plots
          (log-normalized expression), grouped as:
            row 1  neuronal:   Snap25, Rbfox3, Slc17a7, Gad1
            row 2  astrocyte:  Slc1a3, Aqp4, Gfap, Aldh1l1
            row 3  MCT/acetate axis: Slc16a1 (MCT1), Slc16a3, Acss1, Acss2
          The payoff: the acetate-handling machinery (MCT / Acss) localizes to
          the astrocyte island, the premise of the model.

Embedding: PCA(50) on log-normalized expression -> neighbours(15) -> UMAP, all
random_state=0. Cells restricted to the four acute-design groups (H2O_veh,
H2O_MCT1i, EtOH_veh, EtOH_MCT1i) so the map matches the rest of the figures.

NOTE: cell types come from curated marker-gene module scores (master build), not
de-novo clustering, so the UMAP is illustrative of separation rather than an
independent clustering result.

Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso). All paths via
CLI. Vector PDF, editable text (pdf.fonttype=42); point clouds rasterized to keep
the file small while all labels stay as text. One source CSV per version.
"""
from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import anndata as ad
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# ----------------------------- parameters ----------------------------------
GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i']   # acute design

# Display cell-type categories (plot order = back -> front; smaller pops on top)
CT_ORDER = ['Other', 'Oligodendrocyte', 'Microglia', 'Astrocyte',
            'Excitatory', 'Inhibitory']
CT_COLORS = {
    'Excitatory':      '#1f77b4',   # blue
    'Inhibitory':      '#d62728',   # red
    'Astrocyte':       '#2ca02c',   # green
    'Oligodendrocyte': '#9467bd',   # purple
    'Microglia':       '#8c564b',   # brown
    'Other':           '#d9d9d9',   # light gray
}
LEGEND_ORDER = ['Excitatory', 'Inhibitory', 'Astrocyte',
                'Oligodendrocyte', 'Microglia', 'Other']

# Feature-plot grid (3 rows x 4 cols): neuronal / astrocyte / MCT-acetate axis
FEATURE_ROWS = [
    ('Neuronal',          ['Snap25', 'Rbfox3', 'Slc17a7', 'Gad1']),
    ('Astrocyte',         ['Slc1a3', 'Aqp4', 'Gfap', 'Aldh1l1']),
    ('MCT / acetate axis', ['Slc16a1', 'Slc16a7', 'Slc16a3', 'Acss1', 'Acss2']),
]
FEATURE_GENES = [g for _sec, gs in FEATURE_ROWS for g in gs]
ALIAS = {'Slc16a1': 'MCT1', 'Slc16a7': 'MCT2', 'Slc16a3': 'MCT4',
         'Slc1a3': 'GLAST', 'Rbfox3': 'NeuN', 'Slc17a7': 'VGLUT1'}

N_PCS = 50
N_NEIGHBORS = 15
RANDOM_STATE = 0

# ----------------------------- loaders (master build) ------------------------
LABEL_MAP = {
    'Chronic EtOH 1': ('ChronicEtOH', 1), 'MCT1i+EtOH 1': ('EtOH_MCT1i', 1),
    'MCT1i+Veh 1': ('H2O_MCT1i', 1), 'MAT2A CM 1': ('MAT2A_CM', 1),
    'MAT2A CM 2': ('MAT2A_CM', 2), 'MAT2A CM 3': ('MAT2A_CM', 3),
    'MAT2A OE 1': ('MAT2A_OE', 1), 'MAT2A OE 2': ('MAT2A_OE', 2),
    'MAT2A OE 3': ('MAT2A_OE', 3), 'Veh+EtOH 1': ('EtOH_veh', 1),
    'veh +EtOH 2': ('EtOH_veh', 2), 'veh+H20 1': ('H2O_veh', 1),
    'EtOH +MCT1i 2': ('EtOH_MCT1i', 2), 'EtOH +MCT1i 3': ('EtOH_MCT1i', 3),
    'H20+veh 2': ('H2O_veh', 2), 'H20 +MCT1i 2': ('H2O_MCT1i', 2),
    'Chronic EtOH 2': ('ChronicEtOH', 2), 'chronic EtOH 3': ('ChronicEtOH', 3),
    'EtOH+MCT1i 1': ('EtOH_MCT1i', 1), 'H20+MCT1i 1': ('H2O_MCT1i', 1),
    'H20 +MCT1i 1': ('H2O_MCT1i', 1), 'H20 +veh 1': ('H2O_veh', 1),
    'EtOH +veh 1': ('EtOH_veh', 1), 'APP 1': ('APP', 1),
    'APP 2': ('APP', 2), 'APP 3': ('APP', 3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']
INHIB_MARKERS = ['Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln', 'Lhx6']


def load_polygons(csv):
    df = pd.read_csv(csv, comment='#')
    return {name: sub[['X', 'Y']].values.astype(float)
            for name, sub in df.groupby('Selection', sort=False)}


def assign(cells, polygons):
    out = pd.Series(np.full(len(cells), np.nan, dtype=object), index=cells.index)
    pts = cells[['x_centroid', 'y_centroid']].values
    for name, xy in polygons.items():
        inside = MplPath(xy).contains_points(pts)
        out.values[inside & ~out.notna().values] = name
    return out


def process_slide(slide_dir, ann_csv, label):
    cells = pd.read_csv(slide_dir / 'cells.csv.gz').set_index('cell_id')
    polys = load_polygons(ann_csv)
    cells['roi'] = assign(cells, polys)
    cells['group'] = cells['roi'].map(lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)
    cells['slide'] = label
    cells_in = cells.dropna(subset=['group']).copy()
    adata = sc.read_10x_h5(str(slide_dir / 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    adata.obs_names = [f'{label}::{x}' for x in adata.obs_names]
    cells_in.index = [f'{label}::{x}' for x in cells_in.index]
    common = adata.obs_names.intersection(cells_in.index)
    adata = adata[common].copy()
    adata.obs = cells_in.loc[common, ['roi', 'group', 'replicate', 'sample_id', 'slide']].copy()
    return adata


def cell_type_data(adata):
    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gs_present, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    scores = adata.obs[[f'score_{ct}' for ct in MARKERS]]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_ct[scores.max(axis=1) <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best_ct, categories=list(MARKERS.keys()) + ['Unclassified'])
    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory', random_state=0, n_bins=25)
    inhib_present = [g for g in INHIB_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=inhib_present, score_name='score_Inhibitory', random_state=0, n_bins=25)
    is_neuron = adata.obs['celltype'].astype(str) == 'Neuron'
    is_excit = is_neuron & (adata.obs['score_Excitatory'] > 0)
    is_inhib = is_neuron & (adata.obs['score_Inhibitory'] > 0)
    is_astro = adata.obs['celltype'].astype(str) == 'Astrocyte'
    subtype = pd.Series('Other', index=adata.obs.index)
    subtype[is_inhib] = 'Inhibitory'
    subtype[is_excit] = 'Excitatory'  # excit precedence
    subtype[is_astro] = 'Astrocyte'
    adata.obs['subtype'] = pd.Categorical(subtype, categories=['Excitatory', 'Inhibitory', 'Astrocyte', 'Other'])
    return adata


def display_celltype(obs):
    """Collapse subtype + celltype into one display label for the UMAP."""
    sub = obs['subtype'].astype(str)
    cty = obs['celltype'].astype(str)
    out = pd.Series('Other', index=obs.index)
    out[cty == 'Oligodendrocyte'] = 'Oligodendrocyte'
    out[cty == 'Microglia'] = 'Microglia'
    out[sub == 'Astrocyte'] = 'Astrocyte'
    out[sub == 'Inhibitory'] = 'Inhibitory'
    out[sub == 'Excitatory'] = 'Excitatory'
    return out


def process_version(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    ada_a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    ada_b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)          # adata.X is now log-normalized expression
    adata = cell_type_data(adata)
    # restrict to the four acute-design groups
    adata = adata[adata.obs['group'].astype(str).isin(GROUPS)].copy()
    adata.obs['celltype_display'] = pd.Categorical(
        display_celltype(adata.obs), categories=LEGEND_ORDER)
    return adata


def embed(adata):
    """PCA -> neighbours -> UMAP on log-normalized expression (no scaling, sparse)."""
    sc.pp.pca(adata, n_comps=N_PCS, svd_solver='arpack', random_state=RANDOM_STATE)
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    return adata


def dense_col(adata, gene):
    X = adata[:, gene].X
    return np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()


def disp(g):
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


# ----------------------------- drawing --------------------------------------
def draw_celltype_page(pdf, adata, version):
    um = adata.obsm['X_umap']
    lab = adata.obs['celltype_display'].astype(str).values
    fig, ax = plt.subplots(figsize=(8.0, 6.4))
    counts = pd.Series(lab).value_counts().to_dict()
    for ct in CT_ORDER:
        m = lab == ct
        if not m.any():
            continue
        ax.scatter(um[m, 0], um[m, 1], s=2.2, c=CT_COLORS[ct], linewidths=0,
                   rasterized=True)
    handles = [plt.Line2D([0], [0], marker='o', linestyle='', markersize=7,
                          markerfacecolor=CT_COLORS[ct], markeredgecolor='none',
                          label=f'{ct}  (n={counts.get(ct, 0):,})')
               for ct in LEGEND_ORDER if counts.get(ct, 0) > 0]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=8, handletextpad=0.3, labelspacing=0.7)
    ax.set_xlabel('UMAP 1', fontsize=9)
    ax.set_ylabel('UMAP 2', fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.set_title(f'Slides A+B pooled — cell-type UMAP — {version}\n'
                 f'(n={len(lab):,} cells; four acute-design groups; '
                 f'cell types from curated marker scores)',
                 fontsize=10, fontweight='bold')
    pdf.savefig(fig, bbox_inches='tight', dpi=300)
    plt.close(fig)


def draw_feature_page(pdf, adata, version):
    um = adata.obsm['X_umap']
    var = set(adata.var_names)
    nrow, ncol = len(FEATURE_ROWS), max(len(gs) for _s, gs in FEATURE_ROWS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow))
    for ri, (sec, genes) in enumerate(FEATURE_ROWS):
        for ci in range(ncol):
            ax = axes[ri, ci]
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if ci >= len(genes):
                ax.axis('off'); continue
            gene = genes[ci]
            if gene not in var:
                ax.text(0.5, 0.5, f'{disp(gene)}\n(not in panel)', ha='center',
                        va='center', fontsize=8, color='0.5', transform=ax.transAxes)
                ax.set_title(disp(gene), fontsize=9, color='0.5')
                continue
            expr = dense_col(adata, gene)
            order = np.argsort(expr)            # low first -> high on top
            vmax = np.percentile(expr, 99)
            if not np.isfinite(vmax) or vmax <= 0:
                vmax = float(expr.max()) if expr.max() > 0 else 1.0
            scat = ax.scatter(um[order, 0], um[order, 1], s=1.4,
                              c=expr[order], cmap='viridis', vmin=0, vmax=vmax,
                              linewidths=0, rasterized=True)
            ax.set_title(disp(gene), fontsize=9, fontweight='bold')
            cb = fig.colorbar(scat, ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=6)
        axes[ri, 0].set_ylabel(sec, fontsize=9, fontweight='bold')
    fig.suptitle(f'Slides A+B pooled — marker / MCT-axis feature plots — {version}\n'
                 'colour = log-normalized expression (per-gene scale, 0–99th pct); same UMAP as cell-type page',
                 fontsize=10, y=1.0)
    pdf.savefig(fig, bbox_inches='tight', dpi=300)
    plt.close(fig)


def write_source_csv(adata, path):
    um = adata.obsm['X_umap']
    df = pd.DataFrame({
        'cell_id': adata.obs_names,
        'slide': adata.obs['slide'].astype(str).values,
        'group': adata.obs['group'].astype(str).values,
        'celltype_display': adata.obs['celltype_display'].astype(str).values,
        'UMAP1': um[:, 0], 'UMAP2': um[:, 1],
    })
    for g in FEATURE_GENES:
        if g in set(adata.var_names):
            df[f'expr_{g}'] = dense_col(adata, g)
    df.to_csv(path, index=False)


def cover(title, body, path):
    d = fitz.open(); p = d.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 54, 738, 150), title, fontsize=16,
                     fontname='Helvetica-Bold', align=1)
    p.insert_textbox(fitz.Rect(54, 168, 738, 560), body, fontsize=11,
                     fontname='Helvetica')
    d.save(str(path)); d.close()


# ----------------------------- main -----------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-b-ann', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)s | %(message)s')
    today = date.today().isoformat()
    args.outdir.mkdir(parents=True, exist_ok=True)

    tmp = args.outdir / '_tmp_umap'; tmp.mkdir(exist_ok=True)
    pages_pdf = tmp / 'pages.pdf'

    versions = {
        'V1': (args.v1_slide_a_ann, args.v1_slide_b_ann),
        'V2': (args.v2_slide_a_ann, args.v2_slide_b_ann),
    }
    with PdfPages(pages_pdf) as pdf:
        for version, (a_ann, b_ann) in versions.items():
            logging.info('=== Processing %s (pool A+B) ===', version)
            adata = process_version(args.slide_a_dir, a_ann, args.slide_b_dir, b_ann)
            comp = adata.obs['celltype_display'].astype(str).value_counts().to_dict()
            logging.info('%s cell-type composition: %s', version, comp)
            logging.info('%s embedding %d cells x %d genes ...', version,
                         adata.n_obs, adata.n_vars)
            adata = embed(adata)
            csv_path = args.outdir / f'SlidesAB_UMAP_celltype_markers_{version}_{today}.csv'
            write_source_csv(adata, csv_path)
            logging.info('Wrote %s', csv_path.name)
            draw_celltype_page(pdf, adata, version)
            draw_feature_page(pdf, adata, version)
            logging.info('Pages done: %s', version)

    cover_path = tmp / 'cover.pdf'
    cover('Xenium May 2026 — cell-type UMAP + MCT/acetate feature plots',
          f"Naomi Kassahun · {today} · MEWS Lab\n\n"
          "Neurons-vs-astrocytes UMAP for the acute ethanol ± MCT1 inhibition\n"
          "design. Two pages per version:\n\n"
          "  PAGE A  cell-type UMAP — one embedding coloured by cell type\n"
          "    (Excitatory / Inhibitory / Astrocyte / Oligodendrocyte /\n"
          "    Microglia / Other). Neurons and astrocytes occupy distinct\n"
          "    transcriptional territory.\n"
          "  PAGE B  the SAME embedding as a 3x4 grid of feature plots\n"
          "    (log-normalized expression): neuronal markers (Snap25, Rbfox3,\n"
          "    Slc17a7, Gad1), astrocyte markers (Slc1a3, Aqp4, Gfap, Aldh1l1),\n"
          "    and the MCT / acetate axis (Slc16a1/MCT1, Slc16a3, Acss1, Acss2).\n"
          "    The acetate-handling machinery localizes to the astrocyte island.\n\n"
          "Embedding: PCA(50) on log-normalized expression -> neighbours(15) ->\n"
          "UMAP (all random_state=0). Cells restricted to the four acute-design\n"
          "groups (H2O_veh, H2O_MCT1i, EtOH_veh, EtOH_MCT1i); Slides A+B pooled.\n\n"
          "Built for V1 (broad-ROI punches) and V2 (HPC-restricted lasso).\n\n"
          "Caveat: cell types come from curated marker-gene module scores, not\n"
          "de-novo clustering — the UMAP illustrates separation rather than\n"
          "providing an independent clustering result.",
          cover_path)

    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    src = fitz.open(str(pages_pdf))
    out.insert_pdf(src, from_page=0, to_page=2 * len(versions) - 1)
    out_pdf = args.outdir / f'Xenium_UMAP_celltype_markers_V1V2_Summary_{today}.pdf'
    out.save(str(out_pdf)); out.close(); src.close()
    logging.info('Wrote %s (%d pages)', out_pdf, 1 + 2 * len(versions))

    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    logging.info('Done.')


if __name__ == '__main__':
    main()
