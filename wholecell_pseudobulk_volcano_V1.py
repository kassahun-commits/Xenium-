#!/usr/bin/env python3
"""
WHOLE-CELL (region V1) pseudobulk DESeq2 volcano plots — Neuron & Astrocyte.

Purpose
-------
Re-make the V1 Neuron-vs-Astrocyte volcano analysis using *pseudobulk DESeq2*
instead of the cell-level Wilcoxon used in the original V1 master summary
(build_PI_master_NvsA_lfc0.5.py). The original explicitly noted its adj
p-values are "inflated by pseudoreplication"; this version is the replicate-
level (mouse = experimental unit) companion so the two can be compared.

Pipeline is IDENTICAL to the original V1 build, except for the DE test:
  - same broad-ROI ("full punch") V1 annotations,
  - same WHOLE-CELL counts (cell_feature_matrix.h5 = all transcripts/cell),
  - same QC (min_counts=10, min_cells=5), same marker-based cell typing,
  - same 7 contrasts.
Only difference: sum raw counts per (cell type x mouse) -> DESeq2 ~group.

V1 replicate structure (mice/group): H2O_veh 2, H2O_MCT1i 2, EtOH_veh 2,
EtOH_MCT1i 3, ChronicEtOH 3, MAT2A_CM 3, MAT2A_OE 3, APP 3.

Thresholds: padj < 0.05 (proper pseudobulk cutoff) & |log2FC| > 0.5.
(The original Wilcoxon figure used padj < 1e-3; pseudobulk is better calibrated
so the standard 0.05 is appropriate. Both are reported for comparison.)

Reads all paths from CLI (no hardcoding). Editable-text PDFs (fonttype 42),
source-data CSVs alongside.
"""
import argparse
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# ---- constants copied verbatim from build_PI_master_NvsA_lfc0.5.py ----
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
    'EtOH +veh 1': ('EtOH_veh', 1), 'APP 1': ('APP', 1), 'APP 2': ('APP', 2),
    'APP 3': ('APP', 3),
}
MARKERS = {
    'Neuron':          ['Rbfox3', 'Snap25', 'Syn1', 'Syt1', 'Stmn2', 'Map2', 'Tubb3'],
    'Astrocyte':       ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b', 'Aldoc', 'Gja1'],
    'Oligodendrocyte': ['Mog', 'Olig1', 'Olig2', 'Sox10'],
    'Microglia':       ['Cx3cr1', 'Tmem119', 'Csf1r', 'Aif1', 'Trem2'],
}
SUBTYPES = ['Neuron', 'Astrocyte']
SUBTYPE_COLORS = {'Neuron': '#1f77b4', 'Astrocyte': '#d62728'}

# same 7 contrasts as the V1 master (test, reference, label)
COMPARISONS = [
    ('EtOH_veh',    'H2O_veh',  'EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP vs H2O_veh'),
]

LFC_THRESH = 0.5
PADJ_THRESH = 0.05
MIN_TOTAL_COUNT = 10
MIN_CELLS_PER_SAMPLE = 25
XLIM = 5.0


# ---- data pipeline (mirrors the original; keeps RAW counts) ----
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
    cells['group'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][0] if pd.notna(r) and r in LABEL_MAP else None)
    cells['replicate'] = cells['roi'].map(
        lambda r: LABEL_MAP[r][1] if pd.notna(r) and r in LABEL_MAP else None)
    cells['sample_id'] = cells.apply(
        lambda r: f"{r['group']}_{int(r['replicate'])}" if r['group'] else None, axis=1)
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
    best = adata.obs[score_cols].idxmax(axis=1).str.replace('score_', '', regex=False)
    best[adata.obs[score_cols].max(axis=1) <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(best, categories=list(MARKERS.keys()) + ['Unclassified'])
    return adata


def build_v1(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([a, b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    adata.layers['counts'] = adata.X.copy()          # RAW counts for pseudobulk
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)                                # log-norm for cell typing
    adata = cell_type_data(adata)
    return adata


# ---- pseudobulk + DESeq2 ----
def make_pseudobulk(adata_ct):
    counts = adata_ct.layers['counts']
    samples = adata_ct.obs['sample_id'].astype(str).values
    uniq = sorted(pd.unique(samples))
    mat = np.zeros((len(uniq), adata_ct.n_vars), dtype=np.int64)
    n_cells, groups = [], []
    for i, s in enumerate(uniq):
        m = samples == s
        sub = adata_ct[m].layers['counts']
        sub = np.asarray(sub.todense()) if hasattr(sub, 'todense') else np.asarray(sub)
        mat[i] = np.rint(sub.sum(axis=0)).astype(np.int64)
        n_cells.append(int(m.sum()))
        groups.append(adata_ct.obs.loc[m, 'group'].astype(str).iloc[0])
    counts_df = pd.DataFrame(mat, index=uniq, columns=list(adata_ct.var_names))
    meta_df = pd.DataFrame({'sample_id': uniq, 'group': groups, 'n_cells': n_cells}, index=uniq)
    return counts_df, meta_df


def run_deseq(counts_df, meta_df, test, ref):
    keep = meta_df.index[meta_df['group'].isin([test, ref])]
    meta = meta_df.loc[keep].copy()
    gc = meta['group'].value_counts()
    if gc.get(test, 0) < 2 or gc.get(ref, 0) < 2:
        return None
    counts = counts_df.loc[keep]
    counts = counts.loc[:, counts.sum(axis=0) >= MIN_TOTAL_COUNT]
    test_s, ref_s = test.replace('_', '-'), ref.replace('_', '-')
    meta['group'] = pd.Categorical([g.replace('_', '-') for g in meta['group']],
                                   categories=[ref_s, test_s])
    dds = DeseqDataSet(counts=counts, metadata=meta, design_factors='group',
                       ref_level=['group', ref_s], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['group', test_s, ref_s], quiet=True)
    st.summary()
    res = st.results_df.copy()
    res['gene'] = res.index
    return res.reset_index(drop=True)


# ---- plotting ----
def volcano_ax(ax, res, title, color, padj_cut=PADJ_THRESH):
    if res is None or len(res) == 0:
        ax.text(0.5, 0.5, '(insufficient replicates)', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=11, color=color, fontweight='bold')
        ax.set_xticks([]); ax.set_yticks([]); return 0, 0
    d = res.dropna(subset=['log2FoldChange', 'padj']).copy()
    x = d['log2FoldChange'].values
    y = -np.log10(np.clip(d['padj'].values, 1e-300, 1))
    sig = (np.abs(x) > LFC_THRESH) & (d['padj'].values < padj_cut)
    ax.scatter(x[~sig], y[~sig], s=5, c='lightgray', alpha=0.45, edgecolor='none', rasterized=True)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=13, c='#d62728', alpha=0.8, edgecolor='none')
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=13, c='#1f77b4', alpha=0.8, edgecolor='none')
    top = d[sig].reindex(d[sig]['padj'].sort_values().index).head(12)
    for _, r in top.iterrows():
        if -XLIM <= r['log2FoldChange'] <= XLIM:
            ax.annotate(r['gene'], (r['log2FoldChange'], -np.log10(max(r['padj'], 1e-300))),
                        fontsize=6, fontweight='bold', ha='center', va='bottom')
    ax.axhline(-np.log10(padj_cut), ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(LFC_THRESH, ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(-LFC_THRESH, ls='--', color='black', lw=0.4, alpha=0.4)
    ax.set_xlim(-XLIM, XLIM)
    ax.set_xlabel(f'log2 fold change (clipped ±{XLIM:g})', fontsize=8)
    ax.set_ylabel('-log10 adj p', fontsize=8)
    ax.tick_params(labelsize=7)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    n_up = int((sig & (x > 0)).sum()); n_dn = int((sig & (x < 0)).sum())
    ax.text(0.03, 0.97, f'down: {n_dn}', transform=ax.transAxes, fontsize=9, va='top',
            ha='left', fontweight='bold', color='#1f77b4',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#1f77b4', alpha=0.9))
    ax.text(0.97, 0.97, f'up: {n_up}', transform=ax.transAxes, fontsize=9, va='top',
            ha='right', fontweight='bold', color='#d62728',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#d62728', alpha=0.9))
    ax.set_title(title, fontsize=11, color=color, fontweight='bold')
    return n_up, n_dn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slide-a-dir', required=True, type=Path)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--v1-slide-a-ann', required=True, type=Path)
    ap.add_argument('--v1-slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--wilcoxon-de', type=Path, default=None,
                    help='existing V1 Wilcoxon DE long CSV (for the comparison bar chart)')
    ap.add_argument('--date', default=date.today().isoformat())
    args = ap.parse_args()

    out = args.out_dir; out.mkdir(parents=True, exist_ok=True)
    D = args.date

    print('[build] V1 whole-cell object...')
    adata = build_v1(args.slide_a_dir, args.v1_slide_a_ann,
                     args.slide_b_dir, args.v1_slide_b_ann)
    print(f'[build] {adata.n_obs} cells x {adata.n_vars} genes')
    print('[celltype]', adata.obs['celltype'].value_counts().to_dict())

    all_res, pb_meta = [], []
    from matplotlib.backends.backend_pdf import PdfPages
    panel_pdf = out / f'WholeCell_V1_PseudobulkDE_Volcano_NeuronAstrocyte_AllContrasts_{D}.pdf'
    pp = PdfPages(panel_pdf)

    cache = {}
    for ct in SUBTYPES:
        sub = adata[adata.obs['celltype'].astype(str) == ct].copy()
        counts_df, meta_df = make_pseudobulk(sub)
        for _, r in meta_df.iterrows():
            pb_meta.append({'celltype': ct, **r.to_dict()})
        small = meta_df.index[meta_df['n_cells'] < MIN_CELLS_PER_SAMPLE]
        if len(small):
            print(f'  [warn] {ct}: dropping <{MIN_CELLS_PER_SAMPLE}-cell reps {list(small)}')
            counts_df = counts_df.drop(index=small); meta_df = meta_df.drop(index=small)
        for test, ref, label in COMPARISONS:
            res = run_deseq(counts_df, meta_df, test, ref)
            cache[(ct, test, ref)] = res
            if res is not None:
                r2 = res.copy(); r2['celltype'] = ct
                r2['comparison'] = f'{test}_vs_{ref}'
                all_res.append(r2)
                nsig = ((res['padj'] < PADJ_THRESH) & (res['log2FoldChange'].abs() > LFC_THRESH)).sum()
                print(f'  [ok] {ct} {label}: {nsig} sig')
            else:
                print(f'  [skip] {ct} {label}: <2 reps a side')

    # one volcano page per contrast (Neuron | Astrocyte)
    for test, ref, label in COMPARISONS:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
        for ax, ct in zip(axes, SUBTYPES):
            volcano_ax(ax, cache.get((ct, test, ref)), ct, SUBTYPE_COLORS[ct])
        fig.suptitle(f'V1 whole-cell — {label}  (pseudobulk DESeq2; '
                     f'|log2FC|>{LFC_THRESH}, padj<{PADJ_THRESH})',
                     fontsize=12, y=1.01)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pp.savefig(fig, bbox_inches='tight'); plt.close(fig)
    pp.close()
    print(f'[fig] {panel_pdf.name}')

    # results CSV
    if all_res:
        big = pd.concat(all_res, ignore_index=True)
        big['significant'] = (big['padj'] < PADJ_THRESH) & (big['log2FoldChange'].abs() > LFC_THRESH)
        cols = ['celltype', 'comparison', 'gene', 'baseMean', 'log2FoldChange',
                'lfcSE', 'stat', 'pvalue', 'padj', 'significant']
        big = big[[c for c in cols if c in big.columns]]
        big.to_csv(out / f'WholeCell_V1_PseudobulkDE_results_{D}.csv', index=False)
        print(f'[csv] WholeCell_V1_PseudobulkDE_results_{D}.csv ({len(big)} rows)')
    pd.DataFrame(pb_meta).to_csv(out / f'WholeCell_V1_Pseudobulk_replicates_{D}.csv', index=False)

    # comparison bar chart vs the original Wilcoxon
    if args.wilcoxon_de and args.wilcoxon_de.exists() and all_res:
        wx = pd.read_csv(args.wilcoxon_de)
        wx_sig = wx[(wx['logfc'].abs() >= LFC_THRESH) & (wx['padj'] < 1e-3)]
        wx_counts = wx_sig.groupby(['subtype', 'test', 'reference']).size()
        pb_counts = big[big['significant']].groupby(['celltype', 'comparison']).size()
        labels = [lab for _, _, lab in COMPARISONS]
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
        for ax, ct in zip(axes, SUBTYPES):
            wxv = [int(wx_counts.get((ct, t, r), 0)) for t, r, _ in COMPARISONS]
            pbv = [int(pb_counts.get((ct, f'{t}_vs_{r}'), 0)) for t, r, _ in COMPARISONS]
            xx = np.arange(len(labels)); w = 0.38
            b1 = ax.bar(xx - w / 2, wxv, w, color='#d62728', label='Wilcoxon (padj<1e-3)')
            b2 = ax.bar(xx + w / 2, pbv, w, color='#1f77b4', label='Pseudobulk (padj<0.05)')
            ax.bar_label(b1, fontsize=7); ax.bar_label(b2, fontsize=7)
            ax.set_xticks(xx); ax.set_xticklabels([l.replace(' vs ', '\nvs ') for l in labels],
                                                  fontsize=6.5, rotation=0)
            ax.set_title(ct, fontweight='bold'); ax.set_ylabel('# significant genes', fontsize=9)
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
            if ct == 'Neuron':
                ax.legend(fontsize=8, frameon=False)
        fig.suptitle('V1 whole-cell — significant genes: cell-level Wilcoxon vs pseudobulk',
                     fontsize=13)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        cmp_pdf = out / f'WholeCell_V1_DE_WilcoxonVsPseudobulk_{D}.pdf'
        plt.savefig(cmp_pdf, bbox_inches='tight', dpi=200); plt.close(fig)
        print(f'[fig] {cmp_pdf.name}')

    print('[done]')


if __name__ == '__main__':
    main()
