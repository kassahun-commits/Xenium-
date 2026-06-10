#!/usr/bin/env python3
"""
MAT2A OE vs CM — overexpression/stress overlap analysis (V2 HPC, Slide B only).

Question: are the DE genes driven by MAT2A_OE the same as those driven by
MAT2A_CM (catalytic mutant)?  Overexpressing the catalytically-dead mutant can
itself cause cellular stress, so genes shared by BOTH constructs (each vs the
H2O_veh vehicle control) are most likely an overexpression/stress signature,
while OE-unique genes are candidate effects of MAT2A catalytic activity.

Design notes / caveats
  * MAT2A_OE = 1 ROI (1 mouse).  Cell-level Wilcoxon with n=1 animal is
    pseudoreplication -> EXPLORATORY / descriptive only, not replicated.
  * MAT2A_OE, MAT2A_CM and H2O_veh are all on Slide B -> same-slide, no batch
    confound (this is why only Slide B is used).
  * Cell typing, QC and DE reuse the existing V2 pipeline verbatim
    (slidesAB_volcano_grid_V2.py) so parameters stay identical across the
    project: filter_cells min_counts=10, filter_genes min_cells=5,
    normalize_total 1e4, log1p; marker score_genes argmax -> broad celltype.

For each cell type (Neuron, Astrocyte):
  1. DE: MAT2A_OE vs H2O_veh  and  MAT2A_CM vs H2O_veh  (Wilcoxon).
  2. Overlap of significant genes (|log2FC| >= LFC, padj < PADJ):
     shared / OE-unique / CM-unique counts + Jaccard + gene lists.
  3. Heatmap of the top-N significant genes from EACH contrast (union, ~2N
     rows): z-scored mean log-norm expression across H2O_veh / MAT2A_CM /
     MAT2A_OE, rows clustered, a left strip marks each gene's source list.

All paths come from CLI args.  Figures use editable text (pdf.fonttype 42) and
every figure ships with a source-data CSV.
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from matplotlib_venn import venn2, venn2_circles
from scipy.cluster.hierarchy import linkage, leaves_list

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Reuse the V2 pipeline verbatim (cell typing / QC / DE).
sys.path.insert(0, str(Path(__file__).resolve().parent / 'V2_HPC_only'))
import slidesAB_volcano_grid_V2 as base  # noqa: E402

REF = 'H2O_veh'
CONTRASTS = [('MAT2A_OE', REF), ('MAT2A_CM', REF)]
GROUP_ORDER = ['H2O_veh', 'MAT2A_CM', 'MAT2A_OE']
GROUP_LABEL = {'H2O_veh': 'H2O_veh\n(control)',
               'MAT2A_CM': 'MAT2A_CM',
               'MAT2A_OE': 'MAT2A_OE'}
CELLTYPES = ['Neuron', 'Astrocyte']
SRC_COLOR = {'OE only': '#2C7FB8', 'CM only': '#D95F0E', 'shared': '#444444',
             'forced': '#2CA02C'}


def type_cells(adata):
    """Broad marker-score cell typing, identical to slidesAB_volcano_grid_V2.main."""
    for ct, gs in base.MARKERS.items():
        gp = [g for g in gs if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=gp, score_name=f'score_{ct}',
                          random_state=0, n_bins=25)
    score_cols = [f'score_{ct}' for ct in base.MARKERS]
    scores = adata.obs[score_cols]
    best = scores.idxmax(axis=1).str.replace('score_', '', regex=False)
    best[scores.max(axis=1) <= 0] = 'Unclassified'
    adata.obs['celltype'] = pd.Categorical(
        best.values, categories=list(base.MARKERS.keys()) + ['Unclassified'])
    return adata


def top_sig(de, lfc, padj, n):
    """Top-n significant genes ranked by padj asc, |logfc| desc."""
    if de is None:
        return de
    sig = de[(de['logfc'].abs() >= lfc) & (de['padj'] < padj)].copy()
    sig['absfc'] = sig['logfc'].abs()
    sig = sig.sort_values(['padj', 'absfc'], ascending=[True, False])
    return sig.head(n)


def group_mean_expr(adata_ct, genes):
    """Mean log-norm expression per group (rows=genes, cols=GROUP_ORDER)."""
    cols = {}
    for g in GROUP_ORDER:
        sub = adata_ct[adata_ct.obs['group'].astype(str) == g]
        idx = [adata_ct.var_names.get_loc(x) for x in genes]
        X = sub[:, genes].X
        X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
        cols[g] = X.mean(axis=0)
    return pd.DataFrame(cols, index=genes)[GROUP_ORDER]


def draw_fc_bars(ax, fc, title, lfc, padj, clip=8.0):
    """Grouped horizontal bars: OE-vs-ctrl and CM-vs-ctrl log2FC per gene."""
    fc = fc.sort_values('OE_lfc').reset_index(drop=True)  # asc -> biggest on top
    y = np.arange(len(fc))
    oe = fc['OE_lfc'].clip(-clip, clip)
    cm = fc['CM_lfc'].clip(-clip, clip)
    oe_sig = (fc['OE_lfc'].abs() >= lfc) & (fc['OE_padj'] < padj)
    cm_sig = (fc['CM_lfc'].abs() >= lfc) & (fc['CM_padj'] < padj)
    oe_bars = ax.barh(y + 0.2, oe, height=0.4, color='#2C7FB8',
                      label='MAT2A_OE vs control')
    cm_bars = ax.barh(y - 0.2, cm, height=0.4, color='#D95F0E',
                      label='MAT2A_CM vs control')
    # per-bar alpha: significant = opaque, non-significant = pale
    for bar, sig in zip(oe_bars, oe_sig):
        bar.set_alpha(0.95 if sig else 0.4)
    for bar, sig in zip(cm_bars, cm_sig):
        bar.set_alpha(0.95 if sig else 0.4)
    ax.axvline(0, color='k', lw=0.6)
    for x in (lfc, -lfc):
        ax.axvline(x, ls='--', color='gray', lw=0.5, alpha=0.6)
    ax.set_yticks(y)
    labels = [(g + ' *') if src == 'forced' else g
              for g, src in zip(fc['gene'], fc['source'])]
    ax.set_yticklabels(labels, fontsize=4.5)
    for tick, src in zip(ax.get_yticklabels(), fc['source']):
        if src == 'forced':
            tick.set_color('#2CA02C')
    ax.set_ylim(-0.6, len(fc) - 0.4)
    ax.set_xlim(-clip - 0.5, clip + 0.5)
    ax.set_xlabel('log2 fold change vs H2O_veh', fontsize=8)
    ax.tick_params(axis='x', labelsize=7)
    # annotate bars whose true value was clipped
    for i in range(len(fc)):
        for val, off in ((fc['OE_lfc'][i], 0.2), (fc['CM_lfc'][i], -0.2)):
            if abs(val) > clip:
                ax.text(np.sign(val) * clip, y[i] + off, f' {val:.0f}',
                        fontsize=4, va='center',
                        ha='left' if val > 0 else 'right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, loc='lower right', frameon=False)


def draw_venn(ax, oe_set, cm_set, title):
    """Two-circle Venn of OE-vs-ctrl and CM-vs-ctrl DE gene sets."""
    only_oe = oe_set - cm_set
    only_cm = cm_set - oe_set
    shared = oe_set & cm_set
    v = venn2([oe_set, cm_set],
              set_labels=('MAT2A_OE\nvs control', 'MAT2A_CM\nvs control'),
              set_colors=(SRC_COLOR['OE only'], SRC_COLOR['CM only']),
              alpha=0.55, ax=ax)
    # outline circles even when a region is empty
    venn2_circles([oe_set, cm_set], ax=ax, lw=0.8, color='#555555')
    for sid, n in (('10', len(only_oe)), ('01', len(only_cm)),
                   ('11', len(shared))):
        lbl = v.get_label_by_id(sid)
        if lbl is not None:
            lbl.set_text(str(n))
            lbl.set_fontsize(11)
            lbl.set_fontweight('bold')
    for t in v.set_labels:
        if t is not None:
            t.set_fontsize(8)
    ax.set_title(title, fontsize=10, fontweight='bold')
    return only_oe, shared, only_cm


def draw_heatmap(ax, z, genes, sources, title):
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER], fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=4.5)
    ax.tick_params(length=0)
    # left source strip
    for i, src in enumerate(sources):
        ax.add_patch(plt.Rectangle((-0.62, i - 0.5), 0.18, 1.0,
                                    color=SRC_COLOR[src], lw=0,
                                    clip_on=False))
    ax.set_xlim(-0.65, len(GROUP_ORDER) - 0.5)
    ax.set_title(title, fontsize=10, fontweight='bold')
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlideB_V2_MAT2A')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--lfc-thresh', type=float, default=1.0)
    ap.add_argument('--padj-thresh', type=float, default=1e-3)
    ap.add_argument('--top-n', type=int, default=50,
                    help='top significant genes taken from EACH contrast')
    ap.add_argument('--force-genes', default='',
                    help='comma-separated genes ALWAYS added to the heatmap '
                         'regardless of DE significance (e.g. construct tags '
                         'GFP,Mat2a)')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    D, LFC, PADJ, N = args.date, args.lfc_thresh, args.padj_thresh, args.top_n
    force_genes = [g.strip() for g in args.force_genes.split(',') if g.strip()]
    log = open(args.out_dir / f'{args.label}_overlap_log_{D}.txt', 'w')

    def say(*a):
        msg = ' '.join(str(x) for x in a)
        print(msg); log.write(msg + '\n'); log.flush()

    say('=== MAT2A OE/CM stress-overlap (Slide B only) ===')
    say(f'thresholds: |log2FC|>={LFC}, padj<{PADJ:g}; top_n per contrast={N}')

    # ---- load Slide B, QC, type (same params as V2 pipeline) ----
    adata = base.process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = type_cells(adata)
    keep = adata.obs['group'].astype(str).isin(GROUP_ORDER)
    say('Slide B cells (groups of interest):')
    say(adata.obs.loc[keep].groupby(
        ['group', 'celltype'], observed=True).size().to_string())

    de_rows, overlap_rows, hm_pages, fc_pages = [], [], [], []
    venn_pages, venn_rows = [], []
    for ct in CELLTYPES:
        a = adata[adata.obs['celltype'].astype(str) == ct].copy()
        de = {}
        for test, ref in CONTRASTS:
            df, meta = base.run_de(a, test, ref)
            de[test] = df
            say(f'[{ct}] {test} vs {ref}: '
                + ('OK' if df is not None else 'SKIPPED')
                + f"  (test cells={meta['n_test_cells']}, ref cells={meta['n_ref_cells']}, "
                  f"test mice={meta['n_test_mice']}, ref mice={meta['n_ref_mice']})")
            if df is not None:
                d = df.copy()
                d.insert(0, 'celltype', ct)
                d['test'] = test; d['reference'] = ref
                for k, v in meta.items():
                    d[k] = v
                de_rows.append(d)
        if de['MAT2A_OE'] is None or de['MAT2A_CM'] is None:
            say(f'[{ct}] insufficient cells for one contrast - skipping overlap/heatmap')
            continue

        # ---- overlap of significant gene sets ----
        sig = {t: set(de[t][(de[t]['logfc'].abs() >= LFC)
                            & (de[t]['padj'] < PADJ)]['gene']) for t in de}
        oe, cm = sig['MAT2A_OE'], sig['MAT2A_CM']
        shared = oe & cm
        union = oe | cm
        jac = len(shared) / len(union) if union else float('nan')
        overlap_rows.append({
            'celltype': ct, 'n_OE_sig': len(oe), 'n_CM_sig': len(cm),
            'n_shared': len(shared), 'n_OE_only': len(oe - cm),
            'n_CM_only': len(cm - oe), 'jaccard': round(jac, 3),
            'shared_genes': ';'.join(sorted(shared)),
            'OE_only_genes': ';'.join(sorted(oe - cm)),
            'CM_only_genes': ';'.join(sorted(cm - oe))})
        say(f'[{ct}] OVERLAP: OE_sig={len(oe)}, CM_sig={len(cm)}, '
            f'shared={len(shared)}, OE_only={len(oe-cm)}, CM_only={len(cm-oe)}, '
            f'Jaccard={jac:.3f}')

        # ---- direction-split significant sets (for Venn) ----
        def sig_dir(t, sign):
            d = de[t]
            m = (d['padj'] < PADJ) & ((d['logfc'] >= LFC) if sign > 0
                                      else (d['logfc'] <= -LFC))
            return set(d[m]['gene'])
        dir_sets = {(t, s): sig_dir(t, s)
                    for t in ('MAT2A_OE', 'MAT2A_CM') for s in (1, -1)}
        venn_pages.append((ct, dir_sets))
        for s, lbl in ((1, 'up'), (-1, 'down')):
            o, c = dir_sets[('MAT2A_OE', s)], dir_sets[('MAT2A_CM', s)]
            venn_rows.append({
                'celltype': ct, 'direction': lbl,
                'n_OE': len(o), 'n_CM': len(c),
                'n_shared': len(o & c), 'n_OE_only': len(o - c),
                'n_CM_only': len(c - o),
                'shared_genes': ';'.join(sorted(o & c)),
                'OE_only_genes': ';'.join(sorted(o - c)),
                'CM_only_genes': ';'.join(sorted(c - o))})
            say(f'[{ct}] {lbl.upper():4s}: OE={len(o)}, CM={len(c)}, '
                f'shared={len(o & c)}, OE_only={len(o - c)}, CM_only={len(c - o)}')

        # ---- heatmap gene set: top-N sig from each contrast, union ----
        top_oe = top_sig(de['MAT2A_OE'], LFC, PADJ, N)
        top_cm = top_sig(de['MAT2A_CM'], LFC, PADJ, N)
        genes = list(dict.fromkeys(list(top_oe['gene']) + list(top_cm['gene'])))
        in_oe, in_cm = set(top_oe['gene']), set(top_cm['gene'])
        src = ['shared' if (g in in_oe and g in in_cm)
               else ('OE only' if g in in_oe else 'CM only') for g in genes]
        # forced genes (e.g. construct tags) -> always shown regardless of DE
        for fg in force_genes:
            if fg in a.var_names and fg not in genes:
                genes.append(fg); src.append('forced')
            elif fg not in a.var_names:
                say(f'[{ct}] forced gene {fg!r} absent from panel/after QC - skipped')
        if not genes:
            say(f'[{ct}] no significant or forced genes - skipping heatmap')
            continue

        expr = group_mean_expr(a, genes)
        z = expr.sub(expr.mean(axis=1), axis=0).div(
            expr.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
        order = leaves_list(linkage(z.values, method='average')) \
            if len(genes) > 2 else np.arange(len(genes))
        genes_o = [genes[i] for i in order]
        src_o = [src[i] for i in order]
        z_o = z.iloc[order]

        # source CSV alongside figure (per cell type)
        src_csv = args.out_dir / f'{args.label}_Heatmap_{ct}_top{N}each_source_{D}.csv'
        out = z_o.copy()
        out.columns = [f'z_{c}' for c in out.columns]
        for c in GROUP_ORDER:
            out[f'mean_{c}'] = expr.loc[genes_o, c].values
        out.insert(0, 'source_list', src_o)
        for t, dd in (('OE', de['MAT2A_OE']), ('CM', de['MAT2A_CM'])):
            m = dd.set_index('gene')
            out[f'{t}_log2fc'] = m.reindex(genes_o)['logfc'].values
            out[f'{t}_padj'] = m.reindex(genes_o)['padj'].values
        out.index.name = 'gene'
        out.to_csv(src_csv)
        say(f'[{ct}] wrote {src_csv.name}  ({len(genes_o)} genes)')
        hm_pages.append((ct, z_o, genes_o, src_o, len(in_oe), len(in_cm)))

        # fold-change bar data (same gene set)
        oe_m = de['MAT2A_OE'].set_index('gene')
        cm_m = de['MAT2A_CM'].set_index('gene')
        fc = pd.DataFrame({'gene': genes_o, 'source': src_o})
        fc['OE_lfc'] = oe_m.reindex(genes_o)['logfc'].values
        fc['CM_lfc'] = cm_m.reindex(genes_o)['logfc'].values
        fc['OE_padj'] = oe_m.reindex(genes_o)['padj'].values
        fc['CM_padj'] = cm_m.reindex(genes_o)['padj'].values
        fc_pages.append((ct, fc))

    # ---- write DE + overlap tables ----
    if de_rows:
        de_csv = args.out_dir / f'{args.label}_DE_OEvsH2O_CMvsH2O_{D}.csv'
        pd.concat(de_rows, ignore_index=True).to_csv(de_csv, index=False)
        say('wrote', de_csv.name)
    if overlap_rows:
        ov_csv = args.out_dir / f'{args.label}_DEoverlap_summary_{D}.csv'
        pd.DataFrame(overlap_rows).to_csv(ov_csv, index=False)
        say('wrote', ov_csv.name)

    # ---- heatmap PDF (one page per cell type) ----
    if hm_pages:
        pdf_out = args.out_dir / f'{args.label}_Heatmap_top{N}each_OE_CM_{D}.pdf'
        with PdfPages(pdf_out) as pdf:
            for ct, z_o, genes_o, src_o, n_oe, n_cm in hm_pages:
                h = max(4.0, 0.135 * len(genes_o) + 1.6)
                fig, ax = plt.subplots(figsize=(5.4, h))
                title = (f'{ct}: top {n_oe} OE + top {n_cm} CM DE genes '
                         f'(vs H2O_veh)\nz-scored mean log-norm expression')
                im = draw_heatmap(ax, z_o, genes_o, src_o, title)
                cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cb.set_label('row z-score', fontsize=8)
                cb.ax.tick_params(labelsize=7)
                present_src = set(src_o)
                handles = [Patch(facecolor=SRC_COLOR[k], label=k)
                           for k in ['OE only', 'CM only', 'shared', 'forced']
                           if k in present_src]
                ax.legend(handles=handles, title='top-list source',
                          fontsize=7, title_fontsize=7.5,
                          loc='upper left', bbox_to_anchor=(1.25, 1.0),
                          frameon=False)
                fig.text(0.01, 0.005,
                         'EXPLORATORY: MAT2A_OE = 1 ROI (n=1 mouse); '
                         'cell-level Wilcoxon = pseudoreplication.',
                         fontsize=6, color='#888888')
                plt.tight_layout(rect=[0, 0.02, 1, 1])
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        say('wrote', pdf_out.name)

    # ---- fold-change bar chart (two bars/gene: OE-vs-ctrl, CM-vs-ctrl) ----
    if fc_pages:
        fc_pdf = args.out_dir / f'{args.label}_FoldChangeBars_OE_CM_vs_H2O_{D}.pdf'
        for ct, fc in fc_pages:
            fc.to_csv(args.out_dir
                      / f'{args.label}_FoldChangeBars_{ct}_source_{D}.csv',
                      index=False)
        with PdfPages(fc_pdf) as pdf:
            for ct, fc in fc_pages:
                h = max(4.0, 0.135 * len(fc) + 1.6)
                fig, ax = plt.subplots(figsize=(6.0, h))
                draw_fc_bars(ax, fc,
                             f'{ct}: log2FC vs H2O_veh  (* = forced/construct '
                             f'tag; pale bar = not significant)',
                             LFC, PADJ)
                fig.text(0.01, 0.005,
                         'EXPLORATORY: MAT2A_OE = 1 ROI (n=1 mouse); '
                         'cell-level Wilcoxon = pseudoreplication. '
                         'Bars clipped at |log2FC|=8 (true value annotated).',
                         fontsize=6, color='#888888')
                plt.tight_layout(rect=[0, 0.02, 1, 1])
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        say('wrote', fc_pdf.name)

    # ---- Venn diagrams (OE vs CM DE sets, split UP / DOWN per cell type) ----
    if venn_rows:
        venn_csv = args.out_dir / f'{args.label}_Venn_OE_CM_updown_{D}.csv'
        pd.DataFrame(venn_rows).to_csv(venn_csv, index=False)
        say('wrote', venn_csv.name)
    if venn_pages:
        venn_pdf = args.out_dir / f'{args.label}_Venn_OE_CM_updown_{D}.pdf'
        with PdfPages(venn_pdf) as pdf:
            for ct, dir_sets in venn_pages:
                fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6))
                for ax, (s, lbl) in zip(axes, ((1, 'UP'), (-1, 'DOWN'))):
                    draw_venn(ax, dir_sets[('MAT2A_OE', s)],
                              dir_sets[('MAT2A_CM', s)],
                              f'{lbl}-regulated')
                fig.suptitle(
                    f'{ct}: DE gene overlap, OE vs CM (each vs H2O_veh)\n'
                    f'|log2FC| >= {LFC:g}, padj < {PADJ:g}',
                    fontsize=11, fontweight='bold')
                fig.text(0.01, 0.01,
                         'EXPLORATORY: MAT2A_OE = 1 ROI (n=1 mouse); '
                         'cell-level Wilcoxon = pseudoreplication.',
                         fontsize=6, color='#888888')
                plt.tight_layout(rect=[0, 0.03, 1, 0.92])
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        say('wrote', venn_pdf.name)

    log.close()


if __name__ == '__main__':
    main()
