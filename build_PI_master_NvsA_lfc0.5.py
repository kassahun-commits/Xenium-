#!/usr/bin/env python3
"""
Master PI PDF — V1 + V2 — NEURONS vs ASTROCYTES (combined Excit+Inhib) — LFC 0.5.

Differences vs prior master script:
  * Cell types collapsed to two: Neuron (all cells with celltype=='Neuron'
    regardless of excitatory/inhibitory sub-classification) and Astrocyte.
  * LFC significance threshold = 0.5 (instead of 1.0). Volcano dashed lines
    are at +/-0.5, top-hit selection requires |log2FC| >= 0.5.

Layout per contrast (per V1 / V2 section):
  - 1 volcano page (2 cell-type panels: Neuron | Astrocyte)
  - 2 top-hits pages (page 1 = top 25 UP, page 2 = top 25 DOWN) in 5x5 grid
  - Each top-hit panel has all 14 bars (7 condition clusters x Neuron/Astro)

This script runs the full pipeline end-to-end and writes a single combined PDF.
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

ALL_GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i',
              'ChronicEtOH', 'MAT2A_CM', 'MAT2A_OE', 'APP']
SUBTYPES = ['Neuron', 'Astrocyte']  # only two
SUBTYPE_COLORS = {'Neuron':    '#1f77b4',
                  'Astrocyte': '#d62728'}

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


COMPARISONS = [
    ('EtOH_veh',    'H2O_veh',  'Acute alcohol — EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic alcohol — ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'Drug-only — H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'Alcohol + drug — EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'MCT1i rescue — EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2A overexpression — MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP — APP vs H2O_veh'),
]

LFC_THRESH = 0.5
PADJ_THRESH = 1e-3


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


def cell_type_data(adata):
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
    return adata


def process_version(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    """Process one version (V1 or V2). Returns adata with celltype assigned + normalized counts."""
    ada_a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    ada_b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = cell_type_data(adata)
    return adata


def run_de_2subtypes(adata):
    """Run DE for all comparisons in both Neuron and Astrocyte subsets."""
    rows = []
    cache = {}  # (subtype, test, ref) -> (df, meta)
    for subtype in SUBTYPES:
        sub_all = adata[adata.obs['celltype'].astype(str) == subtype].copy()
        for test, ref, _ in COMPARISONS:
            d = sub_all[sub_all.obs['group'].astype(str).isin([test, ref])].copy()
            if d.shape[0] < 30:
                cache[(subtype, test, ref)] = (None, {'n_test_cells': 0, 'n_ref_cells': 0,
                                                       'n_test_mice': 0, 'n_ref_mice': 0})
                continue
            d.obs['de_group'] = d.obs['group'].astype(str)
            n_t = int((d.obs['de_group'] == test).sum())
            n_r = int((d.obs['de_group'] == ref).sum())
            s_t = sorted(d.obs.loc[d.obs['de_group'] == test, 'sample_id'].unique().tolist())
            s_r = sorted(d.obs.loc[d.obs['de_group'] == ref,  'sample_id'].unique().tolist())
            if min(n_t, n_r) < 10:
                cache[(subtype, test, ref)] = (None, {'n_test_cells': n_t, 'n_ref_cells': n_r,
                                                       'n_test_mice': len(s_t), 'n_ref_mice': len(s_r)})
                continue
            sc.tl.rank_genes_groups(d, 'de_group', method='wilcoxon',
                                    reference=ref, n_genes=d.shape[1])
            r = d.uns['rank_genes_groups']
            df = pd.DataFrame({
                'gene':  [str(g) for g in r['names'][test]],
                'logfc': r['logfoldchanges'][test],
                'pval':  r['pvals'][test],
                'padj':  r['pvals_adj'][test],
                'score': r['scores'][test],
            })
            meta = {'n_test_cells': n_t, 'n_ref_cells': n_r,
                    'n_test_mice': len(s_t), 'n_ref_mice': len(s_r)}
            cache[(subtype, test, ref)] = (df, meta)
            df_out = df.copy()
            df_out['subtype'] = subtype
            df_out['test'] = test; df_out['reference'] = ref
            for k, v in meta.items():
                df_out[k] = v
            rows.append(df_out)
    de_long = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return cache, de_long


def compute_persample_means(adata, gene_list):
    persample = {}
    for subtype in SUBTYPES:
        sub = adata[(adata.obs['celltype'].astype(str) == subtype) &
                     (adata.obs['group'].astype(str).isin(ALL_GROUPS))].copy()
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


def plot_volcano(ax, de_df, title, meta, subtitle_color):
    if de_df is None or len(de_df) == 0:
        ax.text(0.5, 0.5, '(insufficient cells)', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=10, color=subtitle_color)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ('top', 'right'): ax.spines[s].set_visible(False)
        return
    x = de_df['logfc'].values
    y = -np.log10(np.clip(de_df['padj'].values, 1e-300, 1))
    sig = (de_df['logfc'].abs() >= LFC_THRESH) & (de_df['padj'] < PADJ_THRESH)
    ax.scatter(x[~sig], y[~sig], s=4, c='lightgray', alpha=0.45, edgecolor='none')
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=12, c='#d62728', alpha=0.8, edgecolor='none')
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=12, c='#1f77b4', alpha=0.8, edgecolor='none')
    top_up = de_df[sig & (de_df['logfc'] > 0)].nlargest(15, 'logfc')
    top_dn = de_df[sig & (de_df['logfc'] < 0)].nsmallest(15, 'logfc')
    for _, row in pd.concat([top_up, top_dn]).iterrows():
        ax.annotate(row['gene'], (row['logfc'], -np.log10(max(row['padj'], 1e-300))),
                    fontsize=6, alpha=0.85)
    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(LFC_THRESH,  ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(-LFC_THRESH, ls='--', color='black', lw=0.4, alpha=0.4)
    ax.set_xlim(-5, 5)  # clip extreme outliers for readability
    ax.set_xlabel('log2 fold change (x-axis clipped to ±5)', fontsize=8)
    ax.set_ylabel('-log10 adj p-value', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    n_up = int((sig & (x > 0)).sum()); n_dn = int((sig & (x < 0)).sum())
    n_info = (f"test cells={meta['n_test_cells']:,} / mice={meta['n_test_mice']}\n"
              f"ref cells={meta['n_ref_cells']:,} / mice={meta['n_ref_mice']}\n"
              f"up={n_up}, down={n_dn}")
    ax.text(0.02, 0.98, n_info, transform=ax.transAxes, fontsize=6,
            va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='gray', alpha=0.85))
    ax.set_title(title, fontsize=10, color=subtitle_color, fontweight='bold')


def write_volcano_pages(de_cache, pdf, version_label):
    for test, ref, label in COMPARISONS:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for ax, subtype in zip(axes, SUBTYPES):
            df, meta = de_cache.get((subtype, test, ref), (None, {}))
            plot_volcano(ax, df, subtype, meta, SUBTYPE_COLORS[subtype])
        fig.suptitle(f'{version_label} — {label}\n'
                     f'(red = up in test, blue = up in reference; '
                     f'sig: |log2FC| >= {LFC_THRESH:g}, padj < {PADJ_THRESH:g})',
                     fontsize=11, y=1.02)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)


def plot_compact_panel(ax, gene, persample_by_subtype, direction_tag=None):
    """Grouped-cell-type bar panel for Neuron + Astrocyte (14 bars per panel)."""
    group_positions = np.arange(len(ALL_GROUPS))
    offsets = {'Neuron': -0.18, 'Astrocyte': 0.18}
    bar_width = 0.32
    y_max = 0.0
    for gi, grp in enumerate(ALL_GROUPS):
        base_x = group_positions[gi]
        for subtype in SUBTYPES:
            d = persample_by_subtype.get(subtype, pd.DataFrame())
            if len(d) == 0: continue
            d2 = d[(d['gene'] == gene) & (d['group'] == grp)]
            if len(d2) == 0: continue
            means = d2['mean_expr_log'].values
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
    fig, axes = plt.subplots(5, 5, figsize=(13.5, 10.5), squeeze=False)
    for i, (gene, dtag) in enumerate(genes_with_dir[:25]):
        r_, c_ = divmod(i, 5)
        plot_compact_panel(axes[r_, c_], gene, persample, direction_tag=dtag)
    for j in range(len(genes_with_dir), 25):
        r_, c_ = divmod(j, 5)
        axes[r_, c_].axis('off')
    legend_handles = [Patch(facecolor=SUBTYPE_COLORS[s], edgecolor='black',
                            alpha=0.85, label=s) for s in SUBTYPES]
    legend_handles.append(plt.Line2D([0], [0], marker='o', color='w',
                                      markerfacecolor='black', markersize=5,
                                      label='per-mouse mean'))
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=0.7, label='SEM'))
    fig.legend(handles=legend_handles, loc='upper center',
               bbox_to_anchor=(0.5, 1.005), ncol=4, frameon=False, fontsize=8)
    fig.suptitle(title, fontsize=11, y=1.025)
    plt.tight_layout(rect=[0, 0, 1, 0.985])
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def top_hits_two_pages(test, ref, label, de_long, persample, pdf,
                        n_per_page=25):
    contrast_de = de_long[(de_long['test'] == test) & (de_long['reference'] == ref) &
                          (de_long['padj'] < PADJ_THRESH)].copy()
    if contrast_de.empty:
        return
    contrast_de['ct_short'] = contrast_de['subtype'].map({'Neuron': 'N', 'Astrocyte': 'A'})

    def best_per_gene(df):
        rows = []
        for gene, sub in df.groupby('gene'):
            tags = []
            for _, r in sub.iterrows():
                tags.append(f"{r['ct_short']}:{'↑' if r['logfc']>0 else '↓'}({r['logfc']:+.1f})")
            rows.append({
                'gene': gene,
                'max_lfc': float(sub['logfc'].max()),
                'min_lfc': float(sub['logfc'].min()),
                'dir_tag': ' '.join(tags),
            })
        return pd.DataFrame(rows)

    up_df = best_per_gene(contrast_de[contrast_de['logfc'] >= LFC_THRESH])
    up_df = up_df.sort_values('max_lfc', ascending=False).head(n_per_page)
    dn_df = best_per_gene(contrast_de[contrast_de['logfc'] <= -LFC_THRESH])
    dn_df = dn_df.sort_values('min_lfc').head(n_per_page)

    plot_top_hits_page(
        [(r['gene'], r['dir_tag']) for _, r in up_df.iterrows()],
        persample,
        f'{label}\nTop UP genes ({len(up_df)} shown, |log2FC|>={LFC_THRESH:g}, padj<{PADJ_THRESH:g})',
        pdf)
    plot_top_hits_page(
        [(r['gene'], r['dir_tag']) for _, r in dn_df.iterrows()],
        persample,
        f'{label}\nTop DOWN genes ({len(dn_df)} shown, |log2FC|>={LFC_THRESH:g}, padj<{PADJ_THRESH:g})',
        pdf)


def make_cover(title, body, out_path):
    cover = fitz.open()
    page = cover.new_page(width=792, height=612)
    page.insert_textbox(fitz.Rect(54, 60, 738, 130), title,
                        fontsize=17, fontname='Helvetica-Bold', align=1)
    page.insert_textbox(fitz.Rect(54, 160, 738, 560), body,
                        fontsize=11, fontname='Helvetica')
    cover.save(str(out_path)); cover.close()


def section_divider(title, out_path):
    sec = fitz.open()
    p = sec.new_page(width=792, height=612)
    p.insert_textbox(fitz.Rect(54, 240, 738, 350), title,
                     fontsize=22, fontname='Helvetica-Bold', align=1)
    sec.save(str(out_path)); sec.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v2-slide-b-ann', required=True, type=Path)
    ap.add_argument('--n-per-page', type=int, default=25)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s')
    logging.info('=== PI master (Neuron vs Astrocyte; LFC %g) ===', LFC_THRESH)

    # ---- V1 ----
    logging.info('Processing V1...')
    v1_adata = process_version(args.slide_a_dir, args.v1_slide_a_ann,
                                args.slide_b_dir, args.v1_slide_b_ann)
    logging.info('V1 cells: %d', v1_adata.shape[0])
    logging.info('Running V1 DE...')
    v1_cache, v1_de = run_de_2subtypes(v1_adata)

    # ---- V2 ----
    logging.info('Processing V2...')
    v2_adata = process_version(args.slide_a_dir, args.v2_slide_a_ann,
                                args.slide_b_dir, args.v2_slide_b_ann)
    logging.info('V2 cells: %d', v2_adata.shape[0])
    logging.info('Running V2 DE...')
    v2_cache, v2_de = run_de_2subtypes(v2_adata)

    # Determine genes to compute means for
    needed = set()
    for de_data in (v1_de, v2_de):
        sub = de_data[(de_data['logfc'].abs() >= LFC_THRESH) & (de_data['padj'] < PADJ_THRESH)]
        for (_test, _ref), grp in sub.groupby(['test', 'reference']):
            up = grp[grp['logfc'] >= LFC_THRESH].nlargest(args.n_per_page * 2, 'logfc')
            dn = grp[grp['logfc'] <= -LFC_THRESH].nsmallest(args.n_per_page * 2, 'logfc')
            needed.update(up['gene'].tolist())
            needed.update(dn['gene'].tolist())
    needed_list = sorted(needed)
    logging.info('Unique top-hit genes: %d', len(needed_list))

    logging.info('Computing V1 per-sample means...')
    v1_persample = compute_persample_means(v1_adata, needed_list)
    logging.info('Computing V2 per-sample means...')
    v2_persample = compute_persample_means(v2_adata, needed_list)

    # ---- Build PDFs ----
    tmp = args.out.parent / '_tmp_master_NvsA'
    tmp.mkdir(parents=True, exist_ok=True)

    # V1 volcano
    v1_volcano_pdf = tmp / 'v1_volcanos.pdf'
    with PdfPages(v1_volcano_pdf) as pdf:
        write_volcano_pages(v1_cache, pdf, 'V1')
    # V1 top hits
    v1_tophits_pdf = tmp / 'v1_tophits.pdf'
    with PdfPages(v1_tophits_pdf) as pdf:
        for test, ref, label in COMPARISONS:
            top_hits_two_pages(test, ref, f'V1 — {label}',
                               v1_de, v1_persample, pdf, args.n_per_page)
    # V2 volcano
    v2_volcano_pdf = tmp / 'v2_volcanos.pdf'
    with PdfPages(v2_volcano_pdf) as pdf:
        write_volcano_pages(v2_cache, pdf, 'V2')
    # V2 top hits
    v2_tophits_pdf = tmp / 'v2_tophits.pdf'
    with PdfPages(v2_tophits_pdf) as pdf:
        for test, ref, label in COMPARISONS:
            top_hits_two_pages(test, ref, f'V2 — {label}',
                               v2_de, v2_persample, pdf, args.n_per_page)
    logging.info('Built component PDFs')

    # Cover + section dividers
    cover_path = tmp / 'cover.pdf'
    cover_body = (
        f"Naomi Kassahun · {date.today().isoformat()} · MEWS Lab\n\n"
        f"Cell types: Neuron (all neurons, excit + inhib combined) and Astrocyte.\n"
        f"Significance thresholds: |log2 FC| >= {LFC_THRESH:g} AND adj p < {PADJ_THRESH:g}.\n\n"
        "Structure:\n"
        "  V1 (broad-ROI / full punches): volcano grid (2 panels per contrast)\n"
        "    + per-contrast top hits — 2 pages per contrast (page 1 = top 25 UP,\n"
        "    page 2 = top 25 DOWN).\n"
        "  V2 (HPC-restricted lasso ROIs): same.\n\n"
        "Top-hit panel layout: 7 condition-group clusters along x. Within each\n"
        "cluster, two bars side-by-side colored by cell type:\n"
        "  blue = Neuron, red = Astrocyte.\n"
        "Bar = mean of per-mouse means; error bar = SEM; dot = individual mouse.\n\n"
        "Direction tag in panel title shows which cell type was the top hit and\n"
        "the log2FC magnitude, e.g. [N:↑(+0.8) A:↓(-0.6)]\n"
        "  N = Neuron, A = Astrocyte.\n\n"
        "V1 samples: H2O_veh=2, H2O_MCT1i=2, EtOH_veh=2, EtOH_MCT1i=3,\n"
        "ChronicEtOH=3, MAT2A_CM=3, MAT2A_OE=3, APP=3 (21 mice).\n"
        "V2 samples: H2O_veh=1, H2O_MCT1i=1, EtOH_veh=1, EtOH_MCT1i=2,\n"
        "ChronicEtOH=1, MAT2A_CM=3, MAT2A_OE=1, APP=3 (HPC-only lasso ROIs).\n\n"
        "Statistics caveat: cell-level Wilcoxon DE — adj p-values are inflated by\n"
        "pseudoreplication. Direction and log2FC magnitude are trustworthy;\n"
        "treat adj p-values as ranking statistics."
    )
    make_cover('Xenium May 2026 — PI summary  (Neuron vs Astrocyte; LFC 0.5)',
               cover_body, cover_path)

    sec_v1 = tmp / 'sec_v1.pdf'
    section_divider("V1 — broad-ROI analysis\n(full punch outlines)", sec_v1)
    sec_v1_top = tmp / 'sec_v1_top.pdf'
    section_divider("V1 top-hits — bar plots per contrast", sec_v1_top)
    sec_v2 = tmp / 'sec_v2.pdf'
    section_divider("V2 — HPC-restricted lasso ROIs", sec_v2)
    sec_v2_top = tmp / 'sec_v2_top.pdf'
    section_divider("V2 top-hits — bar plots per contrast", sec_v2_top)

    # Merge
    out = fitz.open()
    out.insert_pdf(fitz.open(str(cover_path)))
    out.insert_pdf(fitz.open(str(sec_v1)))
    out.insert_pdf(fitz.open(str(v1_volcano_pdf)))
    out.insert_pdf(fitz.open(str(sec_v1_top)))
    out.insert_pdf(fitz.open(str(v1_tophits_pdf)))
    out.insert_pdf(fitz.open(str(sec_v2)))
    out.insert_pdf(fitz.open(str(v2_volcano_pdf)))
    out.insert_pdf(fitz.open(str(sec_v2_top)))
    out.insert_pdf(fitz.open(str(v2_tophits_pdf)))
    out.save(str(args.out))
    out.close()
    logging.info('Wrote %s (%d pages)', args.out, fitz.open(str(args.out)).page_count)

    # Also save the DE tables for record
    v1_de.to_csv(args.out.parent / f'PI_Master_NvsA_V1_DE_{date.today().isoformat()}.csv', index=False)
    v2_de.to_csv(args.out.parent / f'PI_Master_NvsA_V2_DE_{date.today().isoformat()}.csv', index=False)

    # Clean tmp
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
