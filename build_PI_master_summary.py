#!/usr/bin/env python3
"""
Single master PDF for the PI — combines V1 and V2 (HPC-restricted) analyses.

Structure:
  Cover
  V1 section divider
    V1 volcano grid (1 page per contrast, 7 contrasts)
    Per V1 contrast: 2 condensed pages of top-hits bar plots (page 1 = top UP,
                     page 2 = top DOWN). 25 panels per page in 5x5 grid.
  V2 section divider
    V2 volcano grid (1 page per contrast, 7 contrasts)
    Per V2 contrast: 2 condensed pages of top-hits bar plots

Each top-hits panel: grouped-cell-type bar plot — 7 group clusters x 3 cell
types per cluster (blue=Excit, purple=Inhib, red=Astro), per-mouse means with
SEM, individual mouse dots.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from datetime import date
from pathlib import Path

import anndata as ad
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ----- Combined LABEL_MAP (V1 + V2) -----
LABEL_MAP = {
    'Chronic EtOH 1': ('ChronicEtOH', 1),
    'MCT1i+EtOH 1':   ('EtOH_MCT1i', 1),
    'MCT1i+Veh 1':    ('H2O_MCT1i',  1),
    'MAT2A CM 1':     ('MAT2A_CM',   1),
    'MAT2A CM 2':     ('MAT2A_CM',   2),
    'MAT2A CM 3':     ('MAT2A_CM',   3),
    'MAT2A OE 1':     ('MAT2A_OE',   1),
    'MAT2A OE 2':     ('MAT2A_OE',   2),
    'MAT2A OE 3':     ('MAT2A_OE',   3),
    'Veh+EtOH 1':     ('EtOH_veh',   1),
    'veh +EtOH 2':    ('EtOH_veh',   2),
    'veh+H20 1':      ('H2O_veh',    1),
    'EtOH +MCT1i 2':  ('EtOH_MCT1i', 2),
    'EtOH +MCT1i 3':  ('EtOH_MCT1i', 3),
    'H20+veh 2':      ('H2O_veh',    2),
    'H20 +MCT1i 2':   ('H2O_MCT1i',  2),
    'Chronic EtOH 2': ('ChronicEtOH', 2),
    'chronic EtOH 3': ('ChronicEtOH', 3),
    'EtOH+MCT1i 1':   ('EtOH_MCT1i', 1),
    'H20+MCT1i 1':    ('H2O_MCT1i',  1),
    'H20 +MCT1i 1':   ('H2O_MCT1i',  1),
    'H20 +veh 1':     ('H2O_veh',    1),
    'EtOH +veh 1':    ('EtOH_veh',   1),
    'APP 1':          ('APP',        1),
    'APP 2':          ('APP',        2),
    'APP 3':          ('APP',        3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
EXCIT_MARKERS = ['Slc17a7', 'Camk2a', 'Camk2b', 'Satb2', 'Tbr1', 'Neurod6']
INHIB_MARKERS = ['Gad1', 'Gad2', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'Reln', 'Lhx6']

ALL_GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
              'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE', 'APP']
SUBTYPES = ['Excitatory', 'Inhibitory', 'Astrocyte']
SUBTYPE_COLORS = {'Excitatory': '#1f77b4',
                  'Inhibitory': '#9467bd',
                  'Astrocyte':  '#d62728'}

ALIAS = {
    'Slc16a1': 'MCT1', 'Slc16a7': 'MCT2', 'Slc16a3': 'MCT4',
    'Slc17a7': 'VGLUT1', 'Slc17a6': 'VGLUT2',
    'Slc2a1': 'GLUT1', 'Slc2a2': 'GLUT2', 'Slc2a3': 'GLUT3',
    'Slc2a4': 'GLUT4', 'Slc2a5': 'GLUT5', 'Slc5a2': 'SGLT2',
    'Slc1a3': 'GLAST', 'Slc1a1': 'EAAC1', 'Slc1a2': 'GLT-1',
    'Slc32a1': 'VGAT', 'Slc6a11': 'GAT-3', 'Rbfox3': 'NeuN',
    'Crebbp': 'CBP', 'Ep300': 'p300',
    'Kat2a': 'GCN5', 'Kat2b': 'PCAF', 'Kat5': 'Tip60',
    'GFP': 'reporter',
    'Cat': 'catalase', 'Glul': 'glutamine synthetase',
    'Nr3c1': 'GR', 'Nr3c2': 'MR',
    'Atf4': 'ATF4', 'Atf6': 'ATF6', 'Ddit3': 'CHOP', 'Hspa5': 'BiP',
    'Eif2ak3': 'PERK', 'Eif2ak2': 'PKR', 'Eif2ak4': 'GCN2',
    'Tsc1': 'TSC1', 'Tsc2': 'TSC2', 'Mtor': 'mTOR',
    'Ntrk1': 'TrkA', 'Ntrk2': 'TrkB', 'Ntrk3': 'TrkC',
    'Ngfr': 'p75NTR', 'Ntf3': 'NT-3', 'Ntf5': 'NT-4/5',
    'Nr4a1': 'Nur77', 'Nr4a2': 'Nurr1', 'Nr4a3': 'NOR1',
}


def disp(g: str) -> str:
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


COMPARISONS_V1 = [
    ('EtOH_veh',    'H2O_veh',  'Acute alcohol — EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic alcohol — ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'Drug-only — H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'Alcohol + drug — EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'MCT1i rescue — EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2A overexpression — MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP — APP vs H2O_veh'),
]
COMPARISONS_V2 = COMPARISONS_V1  # same set (V2 also has APP and MAT2A_OE)


def load_polygons(csv):
    df = pd.read_csv(csv, comment='#')
    out = {}
    for name, sub in df.groupby('Selection', sort=False):
        out[name] = sub[['X', 'Y']].values.astype(float)
    return out


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


def compute_persample_means(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann, gene_list):
    """Run V1 or V2 pipeline + return persample dict keyed by subtype."""
    ada_a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    ada_b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for ct, gs in MARKERS.items():
        gs_present = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gs_present, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    score_cols = [f'score_{ct}' for ct in MARKERS]
    scores = adata.obs[score_cols]
    best_ct = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best_val = scores.max(axis=1)
    best_ct[best_val <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best_ct,
        categories=list(MARKERS.keys()) + ['Unclassified'])
    excit_present = [g for g in EXCIT_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=excit_present, score_name='score_Excitatory',
                      random_state=0, n_bins=25)
    inhib_present = [g for g in INHIB_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=inhib_present, score_name='score_Inhibitory',
                      random_state=0, n_bins=25)
    is_excit = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Excitatory'] > 0))
    is_inhib = ((adata.obs['celltype'].astype(str) == 'Neuron') &
                (adata.obs['score_Inhibitory'] > 0))
    is_astro = (adata.obs['celltype'].astype(str) == 'Astrocyte')

    persample = {}
    for subtype, mask in [('Excitatory', is_excit),
                          ('Inhibitory', is_inhib),
                          ('Astrocyte',  is_astro)]:
        sub = adata[mask].copy()
        sub = sub[sub.obs['group'].astype(str).isin(ALL_GROUPS)].copy()
        sub.obs['group'] = pd.Categorical(sub.obs['group'].astype(str),
                                          categories=ALL_GROUPS, ordered=True)
        g_arr = sub.obs['group'].astype(str).values
        s_arr = sub.obs['sample_id'].astype(str).values
        present = [g for g in gene_list if g in sub.var_names]
        rows = []
        for gene in present:
            X = sub[:, gene].X
            v = np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()
            df = pd.DataFrame({'expr': v, 'group': g_arr, 'sample_id': s_arr})
            mean = df.groupby(['group', 'sample_id'])['expr'].mean().reset_index()
            mean.rename(columns={'expr': 'mean_expr_log'}, inplace=True)
            mean['gene'] = gene
            rows.append(mean)
        persample[subtype] = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return persample


def plot_compact_panel(ax, gene, persample_by_subtype, direction_tag=None):
    """Compact grouped-cell-type bar panel for one gene."""
    group_positions = np.arange(len(ALL_GROUPS))
    offsets = {'Excitatory': -0.27, 'Inhibitory': 0.0, 'Astrocyte': 0.27}
    bar_width = 0.25
    y_max = 0.0
    for gi, grp in enumerate(ALL_GROUPS):
        base_x = group_positions[gi]
        for subtype in SUBTYPES:
            d = persample_by_subtype.get(subtype, pd.DataFrame())
            d = d[(d['gene'] == gene) & (d['group'] == grp)] if len(d) else d
            if len(d) == 0: continue
            means = d['mean_expr_log'].values
            mh = means.mean()
            sem = means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1 else 0.0
            x = base_x + offsets[subtype]
            ax.bar(x, mh, width=bar_width, color=SUBTYPE_COLORS[subtype],
                   edgecolor='black', linewidth=0.3, alpha=0.85, zorder=2)
            if sem > 0:
                ax.errorbar(x, mh, yerr=sem, color='black', capsize=1.0, lw=0.4, zorder=3)
            np.random.seed(hash((gene, subtype, grp)) % (2**32))
            jit = np.random.uniform(-0.05, 0.05, size=len(means))
            ax.scatter([x] * len(means) + jit, means,
                       c='black', s=3, zorder=4, edgecolor='white', linewidth=0.2)
            y_max = max(y_max, mh + sem, float(means.max()) if len(means) else 0)
    for sep in group_positions[:-1] + 0.5:
        ax.axvline(sep, color='lightgray', lw=0.4, alpha=0.6, zorder=1)
    ax.set_xticks(group_positions)
    ax.set_xticklabels(ALL_GROUPS, rotation=45, ha='right', fontsize=4.5)
    ax.tick_params(axis='y', labelsize=4.5)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.set_xlim(-0.55, len(ALL_GROUPS) - 0.45)
    ax.set_ylim(0, max(y_max * 1.20, 0.05))
    title = disp(gene)
    if direction_tag:
        title += f'  [{direction_tag}]'
    ax.set_title(title, fontsize=7, fontweight='bold')


def plot_top_hits_page(genes_with_dir, persample, title, pdf):
    """5x5 grid (25 panels) of bar plots; landscape letter."""
    fig, axes = plt.subplots(5, 5, figsize=(13.5, 10.5), squeeze=False)
    for i, (gene, dtag) in enumerate(genes_with_dir[:25]):
        r_, c_ = divmod(i, 5)
        plot_compact_panel(axes[r_, c_], gene, persample, direction_tag=dtag)
    for j in range(len(genes_with_dir), 25):
        r_, c_ = divmod(j, 5)
        axes[r_, c_].axis('off')
    # legend bottom-right corner of last panel area
    legend_handles = [Patch(facecolor=SUBTYPE_COLORS[s], edgecolor='black',
                            alpha=0.85, label=s) for s in SUBTYPES]
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                      markerfacecolor='black', markersize=5,
                                      label='per-mouse mean'))
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=0.7, label='SEM'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=5, frameon=False, fontsize=8)
    fig.suptitle(title, fontsize=11, y=1.025)
    plt.tight_layout(rect=[0, 0, 1, 0.985])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def top_hits_two_pages(test, ref, label, de_long, persample, pdf,
                        padj_thresh=1e-3, n_per_page=25):
    """For each contrast, write 2 pages: page 1 = top UP, page 2 = top DOWN."""
    contrast_de = de_long[(de_long['test'] == test) & (de_long['reference'] == ref) &
                          (de_long['padj'] < padj_thresh)].copy()
    if contrast_de.empty:
        return
    # Get max/min logfc per gene across cell types + which subtype achieved it
    contrast_de['ct_short'] = contrast_de['subtype'].map(
        {'Excitatory': 'E', 'Inhibitory': 'I', 'Astrocyte': 'A'})

    # Per gene: which cell-type subtypes (with direction) have it as a top hit
    def best_logfc_per_gene(df):
        # Return one row per gene with max_lfc, min_lfc, dir_tag
        rows = []
        for gene, sub in df.groupby('gene'):
            tags = []
            for _, r in sub.iterrows():
                tags.append(f"{r['ct_short']}:{'↑' if r['logfc']>0 else '↓'}({r['logfc']:+.1f})")
            rows.append({
                'gene': gene,
                'max_lfc': float(sub['logfc'].max()),
                'min_lfc': float(sub['logfc'].min()),
                'dir_tag': ' '.join(tags[:3]),  # up to 3 cell-types in tag
            })
        return pd.DataFrame(rows)

    # Top UP: require max_lfc >= 1
    up_df = best_logfc_per_gene(contrast_de[contrast_de['logfc'] > 1])
    up_df = up_df.sort_values('max_lfc', ascending=False).head(n_per_page)
    # Top DOWN: require min_lfc <= -1
    dn_df = best_logfc_per_gene(contrast_de[contrast_de['logfc'] < -1])
    dn_df = dn_df.sort_values('min_lfc').head(n_per_page)

    plot_top_hits_page(
        [(r['gene'], r['dir_tag']) for _, r in up_df.iterrows()],
        persample,
        f'{label}\nTop UP genes  ({len(up_df)} shown, padj < {padj_thresh:g})',
        pdf)
    plot_top_hits_page(
        [(r['gene'], r['dir_tag']) for _, r in dn_df.iterrows()],
        persample,
        f'{label}\nTop DOWN genes  ({len(dn_df)} shown, padj < {padj_thresh:g})',
        pdf)


def make_cover_pdf(title, subtitle, body, out_path):
    cover = fitz.open()
    page = cover.new_page(width=792, height=612)  # landscape letter
    page.insert_textbox(fitz.Rect(54, 60, 738, 130), title,
                        fontsize=18, fontname='Helvetica-Bold', align=1)
    page.insert_textbox(fitz.Rect(54, 135, 738, 165), subtitle,
                        fontsize=12, fontname='Helvetica-Oblique', align=1)
    page.insert_textbox(fitz.Rect(54, 200, 738, 580), body,
                        fontsize=11, fontname='Helvetica')
    cover.save(str(out_path)); cover.close()


def section_divider(title, out_path):
    sec = fitz.open()
    p = sec.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 240, 738, 340), title,
                     fontsize=24, fontname='Helvetica-Bold', align=1)
    sec.save(str(out_path)); sec.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--v1-de-csv', required=True, type=Path)
    ap.add_argument('--v1-volcano-pdf', required=True, type=Path)
    ap.add_argument('--v2-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-b-ann', required=True, type=Path)
    ap.add_argument('--v2-de-csv', required=True, type=Path)
    ap.add_argument('--v2-volcano-pdf', required=True, type=Path)
    ap.add_argument('--padj-thresh', type=float, default=1e-3)
    ap.add_argument('--n-per-page', type=int, default=25)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s')
    logging.info('=== PI master summary ===')

    # Load DE
    v1_de = pd.read_csv(args.v1_de_csv)
    v2_de = pd.read_csv(args.v2_de_csv)

    # Determine genes needed (union of top-N per direction per contrast across V1+V2)
    needed = set()
    for de_data, COMPS in [(v1_de, COMPARISONS_V1), (v2_de, COMPARISONS_V2)]:
        for test, ref, _ in COMPS:
            sub = de_data[(de_data['test'] == test) & (de_data['reference'] == ref) &
                          (de_data['padj'] < args.padj_thresh)]
            up_genes = sub[sub['logfc'] > 1].nlargest(args.n_per_page * 3, 'logfc')['gene']
            dn_genes = sub[sub['logfc'] < -1].nsmallest(args.n_per_page * 3, 'logfc')['gene']
            needed.update(up_genes.tolist())
            needed.update(dn_genes.tolist())
    needed_list = sorted(needed)
    logging.info('Unique genes to compute means for: %d', len(needed_list))

    # Compute per-sample means for V1 and V2 (separately)
    logging.info('Computing V1 per-sample means...')
    v1_persample = compute_persample_means(
        args.slide_a_dir, args.v1_slide_a_ann,
        args.slide_b_dir, args.v1_slide_b_ann,
        needed_list)
    logging.info('Computing V2 per-sample means...')
    v2_persample = compute_persample_means(
        args.slide_a_dir, args.v2_slide_a_ann,
        args.slide_b_dir, args.v2_slide_b_ann,
        needed_list)

    # Generate V1 and V2 top-hits PDFs (2 pages per contrast)
    tmp_dir = args.out.parent / '_tmp_master_build'
    tmp_dir.mkdir(parents=True, exist_ok=True)

    v1_tophits_pdf = tmp_dir / 'v1_tophits.pdf'
    with PdfPages(v1_tophits_pdf) as pdf:
        for test, ref, label in COMPARISONS_V1:
            top_hits_two_pages(test, ref, f'V1 — {label}',
                               v1_de, v1_persample, pdf,
                               args.padj_thresh, args.n_per_page)
    logging.info('Wrote %s', v1_tophits_pdf)

    v2_tophits_pdf = tmp_dir / 'v2_tophits.pdf'
    with PdfPages(v2_tophits_pdf) as pdf:
        for test, ref, label in COMPARISONS_V2:
            top_hits_two_pages(test, ref, f'V2 — {label}',
                               v2_de, v2_persample, pdf,
                               args.padj_thresh, args.n_per_page)
    logging.info('Wrote %s', v2_tophits_pdf)

    # Build covers + dividers
    cover_path = tmp_dir / 'cover.pdf'
    cover_body = (
        "Xenium May 2026 — PI summary deck (V1 + V2)\n\n"
        f"Naomi Kassahun · {date.today().isoformat()} · MEWS Lab\n\n"
        "Structure:\n"
        "  V1 (broad-ROI / full punches): volcano grid + per-contrast top-hits\n"
        "  V2 (HPC-restricted lasso ROIs): same\n\n"
        "Per contrast: 1 volcano-grid page (3 cell-type panels) + 2 condensed pages\n"
        "of top hits. Top hits = top 25 UP genes (page 1) and top 25 DOWN genes\n"
        "(page 2) by largest |log2FC| across cell types (padj < 1e-3).\n\n"
        "Top-hit panels: grouped-cell-type bar plot. Within each panel:\n"
        "  - 7 condition-group clusters along x (H2O_veh, H2O_MCT1i, EtOH_veh,\n"
        "    EtOH_MCT1i, ChronicEtOH, MAT2A_CM, MAT2A_OE, APP)\n"
        "  - 3 bars per cluster colored by cell type (blue=Excitatory,\n"
        "    purple=Inhibitory, red=Astrocyte)\n"
        "  - Bar = mean of per-mouse means, error bar = SEM, dots = individual mice\n"
        "  - Title tag in brackets shows which cell types had the gene as a top hit\n"
        "    and the direction with magnitude, e.g. [E:↑(+1.5) A:↓(-1.2)]\n\n"
        "Statistics caveat: cell-level Wilcoxon DE — adj p-values are inflated by\n"
        "pseudoreplication (cells within a mouse aren't independent). Direction and\n"
        "log2FC magnitude are trustworthy; treat adj p-values as ranking statistics.\n\n"
        "V1 sample sizes: H2O_veh=2, H2O_MCT1i=2, EtOH_veh=2, EtOH_MCT1i=3,\n"
        "ChronicEtOH=3, MAT2A_CM=3, MAT2A_OE=3, APP=3 (21 mice total).\n\n"
        "V2 sample sizes: H2O_veh=1, H2O_MCT1i=1, EtOH_veh=1, EtOH_MCT1i=2,\n"
        "ChronicEtOH=1, MAT2A_CM=3, MAT2A_OE=1, APP=3 (HPC-only lasso ROIs)."
    )
    make_cover_pdf("Xenium May 2026 — PI summary deck",
                   "V1 (broad-ROI) + V2 (HPC-restricted) — volcano plots & top-hit bar plots",
                   cover_body, cover_path)

    sec_v1 = tmp_dir / 'sec_v1.pdf'
    section_divider("V1 — broad-ROI analysis\n(full punch outlines)", sec_v1)
    sec_v1_top = tmp_dir / 'sec_v1_top.pdf'
    section_divider("V1 top-hits — bar plots per contrast", sec_v1_top)
    sec_v2 = tmp_dir / 'sec_v2.pdf'
    section_divider("V2 — HPC-restricted lasso ROIs", sec_v2)
    sec_v2_top = tmp_dir / 'sec_v2_top.pdf'
    section_divider("V2 top-hits — bar plots per contrast", sec_v2_top)

    # Merge
    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    out.insert_pdf(fitz.open(str(sec_v1)))
    out.insert_pdf(fitz.open(str(args.v1_volcano_pdf)))
    out.insert_pdf(fitz.open(str(sec_v1_top)))
    out.insert_pdf(fitz.open(str(v1_tophits_pdf)))
    out.insert_pdf(fitz.open(str(sec_v2)))
    out.insert_pdf(fitz.open(str(args.v2_volcano_pdf)))
    out.insert_pdf(fitz.open(str(sec_v2_top)))
    out.insert_pdf(fitz.open(str(v2_tophits_pdf)))
    out.save(str(args.out))
    out.close()
    logging.info('Wrote %s', args.out)

    # Clean up
    for f in tmp_dir.iterdir():
        f.unlink()
    tmp_dir.rmdir()


if __name__ == '__main__':
    main()
