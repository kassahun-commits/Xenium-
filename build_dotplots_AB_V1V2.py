#!/usr/bin/env python3
"""
Xenium_May2026 — Slides A+B pooled — alcohol cell-type composition / ratio
figures and alcohol-metabolism-axis dot plots.

Remakes two original Slide-B-only figures, but now:
  * pools Slide A + Slide B,
  * for V1 (broad-ROI punches) and V2 (HPC-restricted lasso),
  * in two cell-type schemes: Neuron vs Astrocyte, and Excitatory vs Inhibitory.

Restricted to the 5 alcohol groups (no MAT2A / APP):
    H2O_veh, H2O_MCT1i, EtOH_veh, EtOH_MCT1i, ChronicEtOH

Figures per (version x scheme):
  1. Cell-type composition (stacked fraction) + ratio bars per sample.
  2. Alcohol-axis dot plot: rows = group, cols = alcohol genes,
     one panel per cell type (dot size = fraction expressing, color = mean expr).

Cell typing follows the master-build method: argmax marker score for
Neuron/Astrocyte/Oligodendrocyte/Microglia, then Excit/Inhib split within
neurons (excitatory takes precedence if both score > 0).

All paths via CLI (MEWS Lab rule). Vector PDFs, editable text. Source CSV per figure.
"""
from __future__ import annotations

import argparse
import logging
from datetime import date
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.path import Path as MplPath

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# --- label map (both slides), identical to master build ---
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

ALCOHOL_GROUPS = ['H2O_veh', 'H2O_MCT1i', 'EtOH_veh', 'EtOH_MCT1i', 'ChronicEtOH']
ALCOHOL_AXIS_GENES = ['Adh5', 'Adh1', 'Cat', 'Cyp2e1', 'Aldh2',
                      'Acss1', 'Acss2', 'Slc2a1', 'Slc2a3']
ALIAS = {'Slc2a1': 'GLUT1', 'Slc2a3': 'GLUT3', 'Cat': 'catalase'}

GROUP_COLORS = {
    'H2O_veh': '#7f7f7f', 'H2O_MCT1i': '#17becf', 'EtOH_veh': '#ff7f0e',
    'EtOH_MCT1i': '#2ca02c', 'ChronicEtOH': '#d62728',
}
CT_COLORS = {'Neuron': '#1f77b4', 'Astrocyte': '#d62728',
             'Oligodendrocyte': '#9467bd', 'Microglia': '#2ca02c',
             'Unclassified': '#cccccc'}
SUBTYPE_COMP_COLORS = {'Excitatory': '#1f77b4', 'Inhibitory': '#9467bd',
                       'Astrocyte': '#d62728', 'Other': '#cccccc'}

SCHEMES = {
    'NeuronVsAstrocyte': {
        'col': 'celltype', 'label': 'Neuron vs Astrocyte',
        'comp_cats': ['Neuron', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Unclassified'],
        'comp_colors': CT_COLORS,
        'ratio_num': 'Neuron', 'ratio_den': 'Astrocyte',
        'panels': ['Neuron', 'Astrocyte'],
    },
    'ExcitVsInhib': {
        'col': 'subtype', 'label': 'Excitatory vs Inhibitory',
        'comp_cats': ['Excitatory', 'Inhibitory', 'Astrocyte', 'Other'],
        'comp_colors': SUBTYPE_COMP_COLORS,
        'ratio_num': 'Excitatory', 'ratio_den': 'Inhibitory',
        'panels': ['Excitatory', 'Inhibitory'],
    },
}


def disp(g):
    return f'{g} ({ALIAS[g]})' if g in ALIAS else g


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


def process_version(slide_a_dir, slide_a_ann, slide_b_dir, slide_b_ann):
    ada_a = process_slide(slide_a_dir, slide_a_ann, 'SlideA')
    ada_b = process_slide(slide_b_dir, slide_b_ann, 'SlideB')
    adata = ad.concat([ada_a, ada_b], axis=0, join='inner', merge='same')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return cell_type_data(adata)


def dense_col(adata, gene):
    X = adata[:, gene].X
    return np.asarray(X.todense()).ravel() if hasattr(X, 'todense') else np.asarray(X).ravel()


# ---------------------------------------------------------------------------
# Figure 1: composition (stacked) + ratio bars
# ---------------------------------------------------------------------------
def make_comp_ratio(adata, version, scheme_key, scheme, outdir, today):
    col = scheme['col']
    obs = adata.obs[adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS)].copy()
    obs['ct'] = obs[col].astype(str)
    counts = obs.groupby(['sample_id', 'ct']).size().unstack(fill_value=0)
    for c in scheme['comp_cats']:
        if c not in counts.columns:
            counts[c] = 0
    counts = counts[scheme['comp_cats']]
    samp = counts.reset_index()
    samp['n_total'] = samp[scheme['comp_cats']].sum(axis=1)
    samp = samp[samp['n_total'] > 0].copy()
    samp['group'] = samp['sample_id'].apply(lambda s: s.rsplit('_', 1)[0])
    samp['replicate'] = samp['sample_id'].apply(lambda s: int(s.rsplit('_', 1)[1]))
    slides = (obs.groupby('sample_id')['slide']
                 .apply(lambda s: ','.join(sorted(set(s)))).rename('slides'))
    samp = samp.merge(slides, on='sample_id')
    samp['group'] = pd.Categorical(samp['group'], categories=ALCOHOL_GROUPS, ordered=True)
    samp = samp.sort_values(['group', 'replicate']).reset_index(drop=True)
    num, den = scheme['ratio_num'], scheme['ratio_den']
    samp['ratio'] = samp[num] / samp[den].clip(lower=1)

    csv_path = outdir / f'SlidesAB_Pooled_CellTypeCompRatio_{scheme_key}_{version}_{today}.csv'
    samp.to_csv(csv_path, index=False)

    samples = samp['sample_id'].tolist()
    x = np.arange(len(samples))
    fig, axes = plt.subplots(1, 2, figsize=(max(9, 0.55 * len(samples) + 5), 4.8))

    # left: stacked composition
    bottom = np.zeros(len(samples))
    for c in scheme['comp_cats']:
        vals = (samp[c] / samp['n_total']).values
        axes[0].bar(x, vals, bottom=bottom, width=0.8,
                    color=scheme['comp_colors'].get(c, '#888'), label=c,
                    edgecolor='white', linewidth=0.4)
        bottom += vals
    axes[0].set_xticks(x); axes[0].set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    axes[0].set_ylabel('Fraction of cells'); axes[0].set_ylim(0, 1)
    axes[0].set_title(f'{version} — cell-type composition per sample')
    axes[0].legend(loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)

    # right: ratio bars
    bar_colors = [GROUP_COLORS.get(g, '#444') for g in samp['group'].astype(str)]
    axes[1].bar(x, samp['ratio'].values, color=bar_colors, edgecolor='black', linewidth=0.5)
    axes[1].set_xticks(x); axes[1].set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
    axes[1].set_ylabel(f'{num} / {den} ratio')
    axes[1].set_title(f'{version} — {num} : {den} ratio')
    axes[1].axhline(1.0, ls='--', color='black', lw=0.6, alpha=0.5)
    handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor='black', label=g)
               for g in ALCOHOL_GROUPS if g in samp['group'].astype(str).tolist()]
    axes[1].legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)

    fig.suptitle(f'Slides A+B pooled — {scheme["label"]} — {version}\n'
                 f'(alcohol groups only; {den} clipped to >=1 in ratio)', fontsize=11, y=1.04)
    plt.tight_layout()
    pdf_path = outdir / f'SlidesAB_Pooled_CellTypeCompRatio_{scheme_key}_{version}_{today}.pdf'
    plt.savefig(pdf_path, bbox_inches='tight'); plt.close()
    logging.info('Wrote %s', pdf_path.name)
    return pdf_path


# ---------------------------------------------------------------------------
# Figure 2: alcohol-axis dot plot (panels = cell types)
# ---------------------------------------------------------------------------
def make_dotplot(adata, version, scheme_key, scheme, outdir, today):
    col = scheme['col']; panels = scheme['panels']
    sub = adata[(adata.obs[col].astype(str).isin(panels)) &
                (adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS))].copy()
    genes_present = [g for g in ALCOHOL_AXIS_GENES if g in sub.var_names]
    rows = []
    for gene in genes_present:
        v = dense_col(sub, gene)
        df = pd.DataFrame({'expr': v, 'group': sub.obs['group'].astype(str).values,
                           'ct': sub.obs[col].astype(str).values})
        agg = (df.groupby(['group', 'ct'])
                 .agg(mean_expr=('expr', 'mean'),
                      frac_expr=('expr', lambda z: float((z > 0).mean())),
                      n=('expr', 'size')).reset_index())
        agg['gene'] = gene
        rows.append(agg)
    dot_df = pd.concat(rows, ignore_index=True)
    csv_path = outdir / f'SlidesAB_Pooled_AlcoholAxisDotplot_{scheme_key}_{version}_{today}.csv'
    dot_df.to_csv(csv_path, index=False)

    vmin = float(dot_df['mean_expr'].min()); vmax = float(dot_df['mean_expr'].max())
    n_g = len(genes_present); n_r = len(ALCOHOL_GROUPS)
    fig, axes = plt.subplots(1, len(panels), figsize=(2.6 + len(panels) * (1.0 + n_g * 0.5),
                                                      1.6 + n_r * 0.34), sharey=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, ct in zip(axes, panels):
        d = dot_df[dot_df['ct'] == ct]
        mean_mat = d.pivot(index='group', columns='gene', values='mean_expr').reindex(index=ALCOHOL_GROUPS, columns=genes_present)
        frac_mat = d.pivot(index='group', columns='gene', values='frac_expr').reindex(index=ALCOHOL_GROUPS, columns=genes_present)
        for i, gi in enumerate(ALCOHOL_GROUPS):
            for j, gj in enumerate(genes_present):
                m = mean_mat.iat[i, j]; f = frac_mat.iat[i, j]
                if pd.isna(m) or pd.isna(f):
                    continue
                ax.scatter(j, i, s=(f * 250) + 5, c=[m], cmap='Reds',
                           vmin=vmin, vmax=vmax, edgecolor='black', linewidth=0.4)
        ax.set_xticks(range(n_g)); ax.set_xticklabels([disp(g) for g in genes_present], rotation=45, ha='right', fontsize=8)
        ax.set_yticks(range(n_r)); ax.set_yticklabels(ALCOHOL_GROUPS, fontsize=9)
        ax.set_title(ct, fontsize=10)
        ax.invert_yaxis(); ax.set_xlim(-0.6, n_g - 0.4); ax.set_ylim(n_r - 0.4, -0.6)

    sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.6, pad=0.02)
    cbar.set_label('Mean log-norm expr', fontsize=8)
    # size legend
    size_handles = [plt.scatter([], [], s=(f * 250) + 5, c='gray', edgecolor='black',
                                 linewidth=0.4, label=f'{int(f*100)}%') for f in (0.25, 0.5, 1.0)]
    axes[-1].legend(handles=size_handles, title='frac. expressing', loc='center left',
                    bbox_to_anchor=(1.18, 0.5), frameon=False, fontsize=7, title_fontsize=7)
    fig.suptitle(f'Slides A+B pooled — alcohol-metabolism axis — {scheme["label"]} — {version}\n'
                 '(dot size = fraction expressing; alcohol groups only)', fontsize=10, y=1.02)
    pdf_path = outdir / f'SlidesAB_Pooled_AlcoholAxisDotplot_{scheme_key}_{version}_{today}.pdf'
    plt.savefig(pdf_path, bbox_inches='tight'); plt.close()
    logging.info('Wrote %s', pdf_path.name)
    return pdf_path


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

    versions = {
        'V1': (args.v1_slide_a_ann, args.v1_slide_b_ann),
        'V2': (args.v2_slide_a_ann, args.v2_slide_b_ann),
    }
    for version, (a_ann, b_ann) in versions.items():
        logging.info('=== Processing %s (pool A+B) ===', version)
        adata = process_version(args.slide_a_dir, a_ann, args.slide_b_dir, b_ann)
        n_alc = int(adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS).sum())
        logging.info('%s cells total=%d, in alcohol groups=%d', version, adata.shape[0], n_alc)
        for scheme_key, scheme in SCHEMES.items():
            comp = (adata.obs[adata.obs['group'].astype(str).isin(ALCOHOL_GROUPS)][scheme['col']]
                    .astype(str).value_counts().to_dict())
            logging.info('  %s | %s composition: %s', version, scheme_key, comp)
            make_comp_ratio(adata, version, scheme_key, scheme, args.outdir, today)
            make_dotplot(adata, version, scheme_key, scheme, args.outdir, today)

    logging.info('Done.')


if __name__ == '__main__':
    main()
