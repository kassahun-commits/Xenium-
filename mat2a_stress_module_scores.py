#!/usr/bin/env python3
"""
MAT2A OE/CM — canonical stress-program module scores (V2 HPC, Slide B only).

Question: is overexpression itself stressing the cells?  The cleanest test is
to score canonical cellular-stress gene programs per cell and compare the three
groups (H2O_veh control, MAT2A_CM catalytic mutant, MAT2A_OE).  If BOTH
constructs raise a stress-module score above control, that is direct evidence
that the act of overexpression — not MAT2A catalytic activity — is stressing
the cells (the CM is catalytically dead, so anything it shares with OE is an
overexpression/stress burden).

Modules (curated canonical markers, auto-filtered to genes present in the
panel; coverage is logged):
  * HeatShock_proteostress  - HSP / chaperone (proteotoxic stress)
  * ISR_ATF4                - integrated stress response / ATF4 targets
  * UPR_ERstress            - unfolded-protein response / ER stress
  * Oxidative_NRF2          - NRF2 oxidative-stress response
  * ReactiveAstro           - pan-/A1 reactive-astrocyte inflammation
  * p53_apoptoticStress     - p53 / apoptotic-stress program

Method: sc.tl.score_genes per module on each cell-type subset, then compare
group distributions.  Significance = Mann-Whitney U of each construct vs
H2O_veh (per module, per cell type).

Caveats
  * MAT2A_OE = 1 ROI (1 mouse) -> cell-level stats are pseudoreplication ->
    EXPLORATORY / descriptive only.
  * All three groups are on Slide B (same-slide, no batch confound).
  * QC / cell typing reuse the V2 pipeline verbatim (identical parameters).

All paths come from CLI args.  Figures use editable text (pdf.fonttype 42) and
ship with source-data CSVs.
"""
from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Reuse the overlap script (and through it the V2 pipeline) for identical
# slide processing / QC / cell typing.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mat2a_oe_cm_stress_overlap as ov  # noqa: E402
base = ov.base

GROUP_ORDER = ov.GROUP_ORDER          # ['H2O_veh', 'MAT2A_CM', 'MAT2A_OE']
GROUP_LABEL = ov.GROUP_LABEL
CELLTYPES = ov.CELLTYPES               # ['Neuron', 'Astrocyte']
REF = ov.REF                           # 'H2O_veh'
GROUP_COLOR = {'H2O_veh': '#9E9E9E', 'MAT2A_CM': '#D95F0E', 'MAT2A_OE': '#2C7FB8'}

# Curated canonical stress programs (mouse symbols).  Auto-filtered to panel.
MODULES = {
    'HeatShock_proteostress': ['Hspa1a', 'Hspa1b', 'Hspa5', 'Hsp90aa1',
                               'Hsp90ab1', 'Hspb1', 'Dnajb1', 'Dnaja1',
                               'Hsph1', 'Bag3', 'Hspa8', 'Hspa4', 'Hspd1',
                               'Dnajc3'],
    'ISR_ATF4': ['Atf4', 'Ddit3', 'Trib3', 'Atf3', 'Ppp1r15a', 'Asns',
                 'Chac1', 'Eif2ak3', 'Nupr1', 'Vegfa', 'Slc7a11', 'Sesn2'],
    'UPR_ERstress': ['Xbp1', 'Atf6', 'Ern1', 'Hspa5', 'Herpud1', 'Edem1',
                     'Pdia3', 'Calr', 'Manf', 'Sel1l', 'Dnajb9'],
    'Oxidative_NRF2': ['Nfe2l2', 'Hmox1', 'Nqo1', 'Gclc', 'Gclm', 'Sod1',
                       'Sod2', 'Gpx1', 'Txn1', 'Srxn1', 'Gsr', 'Prdx1'],
    'ReactiveAstro': ['Gfap', 'Vim', 'Serpina3n', 'C4b', 'C3', 'Gbp2',
                      'Cxcl10', 'Lcn2', 'Steap4', 'S100b', 'Aspg', 'Ggta1',
                      'Hspb1'],
    'p53_apoptoticStress': ['Trp53', 'Cdkn1a', 'Bax', 'Bbc3', 'Mdm2',
                            'Gadd45a', 'Phlda3', 'Pmaip1', 'Ccng1', 'Eda2r'],
}


def stars(p):
    if p is None or np.isnan(p):
        return 'n/a'
    return ('***' if p < 1e-3 else '**' if p < 1e-2 else
            '*' if p < 5e-2 else 'ns')


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--slide-b-dir', required=True, type=Path)
    ap.add_argument('--slide-b-ann', required=True, type=Path)
    ap.add_argument('--out-dir', required=True, type=Path)
    ap.add_argument('--label', default='SlideB_V2_MAT2A')
    ap.add_argument('--date', default=date.today().isoformat())
    ap.add_argument('--min-module-genes', type=int, default=3,
                    help='skip a module if fewer than this many genes are '
                         'present in the panel after QC')
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    D = args.date
    log = open(args.out_dir / f'{args.label}_StressModules_log_{D}.txt', 'w')

    def say(*a):
        msg = ' '.join(str(x) for x in a)
        print(msg); log.write(msg + '\n'); log.flush()

    say('=== MAT2A stress-module scoring (Slide B only) ===')

    # ---- load Slide B, QC, type (identical params to V2 pipeline) ----
    adata = base.process_slide(args.slide_b_dir, args.slide_b_ann, 'SlideB')
    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = ov.type_cells(adata)
    keep = adata.obs['group'].astype(str).isin(GROUP_ORDER)
    adata = adata[keep].copy()

    # ---- resolve module gene sets against the panel; log coverage ----
    mod_genes = {}
    for m, gs in MODULES.items():
        present = [g for g in gs if g in adata.var_names]
        say(f'[module] {m}: {len(present)}/{len(gs)} genes present '
            f'-> {present}')
        if len(present) >= args.min_module_genes:
            mod_genes[m] = present
        else:
            say(f'[module] {m}: too few genes (<{args.min_module_genes}) '
                f'- SKIPPED')

    summary_rows, percell_frames = [], []
    page_data = []  # (celltype, {module: {group: array}}, {module: {grp:p}})
    for ct in CELLTYPES:
        a = adata[adata.obs['celltype'].astype(str) == ct].copy()
        grp = a.obs['group'].astype(str).values
        say(f'[{ct}] cells per group: '
            + ', '.join(f'{g}={int((grp == g).sum())}' for g in GROUP_ORDER))
        # score every module on this cell-type subset
        for m, present in mod_genes.items():
            sc.tl.score_genes(a, gene_list=present, score_name=f'mod_{m}',
                              random_state=0, n_bins=25)
        pc = a.obs[['group'] + [f'mod_{m}' for m in mod_genes]].copy()
        pc.columns = ['group'] + list(mod_genes)
        pc.insert(0, 'celltype', ct)
        percell_frames.append(pc.reset_index().rename(
            columns={'index': 'cell_id'}))

        dist, pvals = {}, {}
        for m in mod_genes:
            by = {g: pc.loc[pc['group'] == g, m].values for g in GROUP_ORDER}
            dist[m] = by
            ref = by[REF]
            pvals[m] = {}
            for g in GROUP_ORDER:
                if g == REF:
                    continue
                try:
                    p = mannwhitneyu(by[g], ref, alternative='two-sided').pvalue
                except ValueError:
                    p = np.nan
                pvals[m][g] = p
                summary_rows.append({
                    'celltype': ct, 'module': m, 'group': g,
                    'n_cells': len(by[g]),
                    'mean_score': float(np.mean(by[g])),
                    'sem': float(np.std(by[g], ddof=1) / np.sqrt(len(by[g])))
                    if len(by[g]) > 1 else np.nan,
                    'ref_group': REF,
                    'ref_mean': float(np.mean(ref)),
                    'delta_vs_ref': float(np.mean(by[g]) - np.mean(ref)),
                    'mwu_p_vs_ref': p, 'sig': stars(p)})
                say(f'[{ct}] {m}: {g} mean={np.mean(by[g]):+.4f} '
                    f'(ctrl={np.mean(ref):+.4f}, d={np.mean(by[g])-np.mean(ref):+.4f}) '
                    f'MWU p={p:.2e} {stars(p)}')
            # also record the control group's own mean once
            summary_rows.append({
                'celltype': ct, 'module': m, 'group': REF,
                'n_cells': len(ref), 'mean_score': float(np.mean(ref)),
                'sem': float(np.std(ref, ddof=1) / np.sqrt(len(ref)))
                if len(ref) > 1 else np.nan,
                'ref_group': REF, 'ref_mean': float(np.mean(ref)),
                'delta_vs_ref': 0.0, 'mwu_p_vs_ref': np.nan, 'sig': ''})
        page_data.append((ct, dist, pvals))

    # ---- write source CSVs ----
    if percell_frames:
        pc_csv = args.out_dir / f'{args.label}_StressModules_scores_percell_{D}.csv'
        pd.concat(percell_frames, ignore_index=True).to_csv(pc_csv, index=False)
        say('wrote', pc_csv.name)
    if summary_rows:
        sm_csv = args.out_dir / f'{args.label}_StressModules_summary_{D}.csv'
        pd.DataFrame(summary_rows).to_csv(sm_csv, index=False)
        say('wrote', sm_csv.name)

    # ---- violin figure: one page per cell type, one subplot per module ----
    mods = list(mod_genes)
    n = len(mods)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    pdf_out = args.out_dir / f'{args.label}_StressModules_violins_{D}.pdf'
    with PdfPages(pdf_out) as pdf:
        for ct, dist, pvals in page_data:
            fig, axes = plt.subplots(nrow, ncol,
                                     figsize=(4.0 * ncol, 3.4 * nrow))
            axes = np.atleast_1d(axes).ravel()
            for j, m in enumerate(mods):
                ax = axes[j]
                data = [dist[m][g] for g in GROUP_ORDER]
                parts = ax.violinplot(data, positions=range(len(GROUP_ORDER)),
                                      showmeans=False, showextrema=False,
                                      widths=0.8)
                for body, g in zip(parts['bodies'], GROUP_ORDER):
                    body.set_facecolor(GROUP_COLOR[g])
                    body.set_alpha(0.65)
                    body.set_edgecolor('#333333')
                    body.set_linewidth(0.5)
                # mean +/- sem marker
                for i, g in enumerate(GROUP_ORDER):
                    arr = dist[m][g]
                    mu = np.mean(arr)
                    se = (np.std(arr, ddof=1) / np.sqrt(len(arr))
                          if len(arr) > 1 else 0)
                    ax.errorbar(i, mu, yerr=se, fmt='o', color='black',
                                ms=3.5, lw=1.0, capsize=2.5, zorder=5)
                # significance brackets vs control
                ymax = max(np.percentile(d, 99) for d in data)
                ymin = min(np.percentile(d, 1) for d in data)
                span = (ymax - ymin) or 1.0
                for i, g in enumerate(GROUP_ORDER):
                    if g == REF:
                        continue
                    p = pvals[m].get(g)
                    ax.text(i, ymax + 0.06 * span, stars(p), ha='center',
                            va='bottom', fontsize=8, fontweight='bold')
                ax.set_xticks(range(len(GROUP_ORDER)))
                ax.set_xticklabels([GROUP_LABEL[g].replace('\n', ' ')
                                    for g in GROUP_ORDER], fontsize=7,
                                   rotation=20, ha='right')
                ax.set_title(f'{m}\n({len(mod_genes[m])} genes)', fontsize=8.5,
                             fontweight='bold')
                ax.set_ylabel('module score', fontsize=8)
                ax.tick_params(axis='y', labelsize=7)
                ax.set_ylim(ymin - 0.08 * span, ymax + 0.20 * span)
                for s in ('top', 'right'):
                    ax.spines[s].set_visible(False)
            for j in range(n, len(axes)):
                axes[j].set_visible(False)
            handles = [Patch(facecolor=GROUP_COLOR[g], alpha=0.65,
                             label=GROUP_LABEL[g].replace('\n', ' '))
                       for g in GROUP_ORDER]
            fig.legend(handles=handles, loc='lower center', ncol=3,
                       fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0))
            fig.suptitle(
                f'{ct}: canonical stress-program module scores '
                f'(H2O_veh vs MAT2A_CM vs MAT2A_OE)\n'
                f'* p<0.05  ** p<0.01  *** p<0.001 (Mann-Whitney vs control); '
                f'bars = mean +/- SEM',
                fontsize=10.5, fontweight='bold')
            fig.text(0.01, 0.005,
                     'EXPLORATORY: MAT2A_OE = 1 ROI (n=1 mouse); '
                     'cell-level test = pseudoreplication.',
                     fontsize=6, color='#888888')
            plt.tight_layout(rect=[0, 0.05, 1, 0.93])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    say('wrote', pdf_out.name)
    log.close()


if __name__ == '__main__':
    main()
