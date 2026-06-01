#!/usr/bin/env python3
"""
Splice new volcano-grid pages (with x-axis clipped to ±5) into the existing
master PDFs. Keeps cover, section dividers, and all bar-plot pages unchanged.

Inputs:
  - Existing master PDF (PI_Master_Summary_*)
  - V1 DE CSV (long form: gene, logfc, padj, subtype, test, reference, n_*)
  - V2 DE CSV (same)
  - Cell types (Excit/Inhib/Astro or Neuron/Astro depending on master)
  - LFC threshold (matches the bar-plot threshold of the master)

Output:
  - New master PDF with volcano pages swapped, "_xlim5" suffix added.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Contrast list MUST match what was in the original master's volcano section
COMPARISONS = [
    ('EtOH_veh',    'H2O_veh',  'Acute alcohol — EtOH_veh vs H2O_veh'),
    ('ChronicEtOH', 'H2O_veh',  'Chronic alcohol — ChronicEtOH vs H2O_veh'),
    ('H2O_MCT1i',   'H2O_veh',  'Drug-only — H2O_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'H2O_veh',  'Alcohol + drug — EtOH_MCT1i vs H2O_veh'),
    ('EtOH_MCT1i',  'EtOH_veh', 'MCT1i rescue — EtOH_MCT1i vs EtOH_veh'),
    ('MAT2A_OE',    'MAT2A_CM', 'MAT2A overexpression — MAT2A_OE vs MAT2A_CM'),
    ('APP',         'H2O_veh',  'APP — APP vs H2O_veh'),
]

PADJ_THRESH = 1e-3
XLIM = 5.0

SUBTYPE_COLORS_3 = {'Excitatory': '#1f77b4',
                    'Inhibitory': '#9467bd',
                    'Astrocyte':  '#d62728'}
SUBTYPE_COLORS_2 = {'Neuron':    '#1f77b4',
                    'Astrocyte': '#d62728'}


def plot_volcano(ax, de_df, title, meta, color, lfc_thresh):
    if de_df is None or len(de_df) == 0:
        ax.text(0.5, 0.5, '(insufficient cells)', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=10, color=color)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ('top', 'right'): ax.spines[s].set_visible(False)
        return
    x = de_df['logfc'].values
    y = -np.log10(np.clip(de_df['padj'].values, 1e-300, 1))
    sig = (de_df['logfc'].abs() >= lfc_thresh) & (de_df['padj'] < PADJ_THRESH)
    ax.scatter(x[~sig], y[~sig], s=4, c='lightgray', alpha=0.45, edgecolor='none')
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=12, c='#d62728', alpha=0.8, edgecolor='none')
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=12, c='#1f77b4', alpha=0.8, edgecolor='none')
    top_up = de_df[sig & (de_df['logfc'] > 0)].nlargest(15, 'logfc')
    top_dn = de_df[sig & (de_df['logfc'] < 0)].nsmallest(15, 'logfc')
    for _, row in pd.concat([top_up, top_dn]).iterrows():
        # Only label genes within the xlim range
        lfc = row['logfc']
        if -XLIM <= lfc <= XLIM:
            ax.annotate(row['gene'], (lfc, -np.log10(max(row['padj'], 1e-300))),
                        fontsize=6, alpha=0.85)
    ax.axhline(-np.log10(PADJ_THRESH), ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(lfc_thresh,  ls='--', color='black', lw=0.4, alpha=0.4)
    ax.axvline(-lfc_thresh, ls='--', color='black', lw=0.4, alpha=0.4)
    ax.set_xlim(-XLIM, XLIM)  # <-- clip
    ax.set_xlabel('log2 fold change  (clipped to ±5)', fontsize=8)
    ax.set_ylabel('-log10 adj p-value', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    n_up = int((sig & (x > 0)).sum()); n_dn = int((sig & (x < 0)).sum())
    n_off = int(((de_df['logfc'].abs() > XLIM) & sig).sum())
    n_info = (f"test cells={meta['n_test_cells']:,} / mice={meta['n_test_mice']}\n"
              f"ref cells={meta['n_ref_cells']:,} / mice={meta['n_ref_mice']}\n"
              f"up={n_up}, down={n_dn}"
              + (f"  ({n_off} off-axis)" if n_off > 0 else ""))
    ax.text(0.02, 0.98, n_info, transform=ax.transAxes, fontsize=6,
            va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='gray', alpha=0.85))
    ax.set_title(title, fontsize=10, color=color, fontweight='bold')


def write_volcano_grid_pdf(de_long, subtypes, subtype_colors, version_label,
                            lfc_thresh, out_pdf):
    n_panels = len(subtypes)
    fig_w = 5.0 * n_panels + (0.5 if n_panels == 3 else 1.0)
    with PdfPages(out_pdf) as pdf:
        for test, ref, label in COMPARISONS:
            fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, 5.2))
            if n_panels == 1:
                axes = [axes]
            for ax, subtype in zip(axes, subtypes):
                df = de_long[(de_long['subtype'] == subtype) &
                             (de_long['test'] == test) &
                             (de_long['reference'] == ref)]
                if len(df) == 0:
                    meta = {'n_test_cells': 0, 'n_ref_cells': 0,
                            'n_test_mice': 0, 'n_ref_mice': 0}
                    plot_volcano(ax, None, subtype, meta, subtype_colors[subtype], lfc_thresh)
                else:
                    meta = {
                        'n_test_cells': int(df['n_test_cells'].iloc[0]) if 'n_test_cells' in df.columns else 0,
                        'n_ref_cells':  int(df['n_ref_cells'].iloc[0])  if 'n_ref_cells'  in df.columns else 0,
                        'n_test_mice':  int(df['n_test_mice'].iloc[0])  if 'n_test_mice'  in df.columns else 0,
                        'n_ref_mice':   int(df['n_ref_mice'].iloc[0])   if 'n_ref_mice'   in df.columns else 0,
                    }
                    plot_volcano(ax, df, subtype, meta, subtype_colors[subtype], lfc_thresh)
            fig.suptitle(f'{version_label} — {label}\n'
                         f'(red = up in test, blue = up in reference; '
                         f'sig: |log2FC| >= {lfc_thresh:g}, padj < {PADJ_THRESH:g}; xlim ±{XLIM:g})',
                         fontsize=11, y=1.02)
            plt.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)


def splice_master(orig_pdf_path, new_v1_volcanos, new_v2_volcanos, out_path,
                   v1_volcano_pages=(2, 8), v1_topdiv_to_v2div=(9, 24),
                   v2_volcano_pages=(25, 31), v2_topdiv_to_end=(32, 46)):
    """Reassemble:
      cover + v1_divider (orig pages 0-1)
      + new v1 volcano grid
      + v1 top-hits divider + bar pages (orig pages 9-23)
      + v2 divider (orig page 24)
      + new v2 volcano grid
      + v2 top-hits divider + bar pages (orig pages 32-46)
    """
    src = fitz.open(str(orig_pdf_path))
    nv1 = fitz.open(str(new_v1_volcanos))
    nv2 = fitz.open(str(new_v2_volcanos))
    out = fitz.open()
    # cover + v1 divider
    out.insert_pdf(src, from_page=0, to_page=1)
    # new v1 volcanos
    out.insert_pdf(nv1)
    # v1 top divider + top hits (skip old v1 volcanos)
    out.insert_pdf(src, from_page=v1_topdiv_to_v2div[0], to_page=v1_topdiv_to_v2div[1] - 1)
    # v2 divider
    out.insert_pdf(src, from_page=v1_topdiv_to_v2div[1], to_page=v1_topdiv_to_v2div[1])
    # new v2 volcanos
    out.insert_pdf(nv2)
    # v2 top divider + top hits (skip old v2 volcanos)
    out.insert_pdf(src, from_page=v2_topdiv_to_end[0], to_page=v2_topdiv_to_end[1])
    out.save(str(out_path))
    src.close(); nv1.close(); nv2.close(); out.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--master-in', required=True, type=Path)
    ap.add_argument('--master-out', required=True, type=Path)
    ap.add_argument('--v1-de-csv', required=True, type=Path)
    ap.add_argument('--v2-de-csv', required=True, type=Path)
    ap.add_argument('--lfc-thresh', type=float, required=True)
    ap.add_argument('--subtypes', required=True,
                    help='Comma-separated list of cell types in volcano grid panels; '
                         'either "Excitatory,Inhibitory,Astrocyte" or "Neuron,Astrocyte"')
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s | %(levelname)s | %(message)s')
    subtypes = [s.strip() for s in args.subtypes.split(',')]
    if subtypes == ['Excitatory', 'Inhibitory', 'Astrocyte']:
        colors = SUBTYPE_COLORS_3
    elif subtypes == ['Neuron', 'Astrocyte']:
        colors = SUBTYPE_COLORS_2
    else:
        raise SystemExit(f'Unsupported subtype set: {subtypes}')

    v1_de = pd.read_csv(args.v1_de_csv)
    v2_de = pd.read_csv(args.v2_de_csv)
    logging.info('V1 DE: %d rows; V2 DE: %d rows', len(v1_de), len(v2_de))

    tmp = args.master_out.parent / '_tmp_splice'
    tmp.mkdir(parents=True, exist_ok=True)
    v1_vp = tmp / 'v1_volcanos_xlim5.pdf'
    v2_vp = tmp / 'v2_volcanos_xlim5.pdf'

    write_volcano_grid_pdf(v1_de, subtypes, colors, 'V1', args.lfc_thresh, v1_vp)
    logging.info('Wrote %s', v1_vp)
    write_volcano_grid_pdf(v2_de, subtypes, colors, 'V2', args.lfc_thresh, v2_vp)
    logging.info('Wrote %s', v2_vp)

    splice_master(args.master_in, v1_vp, v2_vp, args.master_out)
    logging.info('Wrote spliced master: %s', args.master_out)

    # Clean tmp
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
