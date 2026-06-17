"""
Candidate Gene Lookup — FC & Significance Across All Fractions & Conditions
============================================================================
Genes of interest (metabolism, alcohol, stress-response, plasticity):
  MCT transporters, GLUTs, PFKs, PDH/PDK, Acss2, Acly,
  Cat, Cyp2e1, Adh5, Aldh2, Grin2b

Pulls pre-computed FC (log2) and corrected p-value from each fraction sheet.
Significance: |FC| > 0.5  AND  corrected p-value > 3.3

Outputs:
  CandidateGenes_table.xlsx
  CandidateGenes_heatmap.pdf / .png
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import matplotlib.ticker as ticker

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_XLSX   = os.path.join(SCRIPT_DIR, 'CandidateGenes_table.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'CandidateGenes_heatmap.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'CandidateGenes_heatmap.png')

# ── Gene list with biological groupings ───────────────────────────────────────
GENE_GROUPS = {
    'MCT\ntransporter':        ['Slc16a1'],
    'Glucose\ntransporters':   ['Slc2a1', 'Slc2a3'],
    'Glycolysis\n(PFK)':       ['Pfkl', 'Pfkm', 'Pfkp'],
    'Pyruvate\nmetabolism':    ['Pdha1', 'Pdk1', 'Pdk2', 'Pdk3'],
    'Acetyl-CoA\nsynthesis':   ['Acss2', 'Acly'],
    'Alcohol\nmetabolism':     ['Adh5', 'Aldh2', 'Cyp2e1'],
    'Oxidative\nstress':       ['Cat'],
    'Neuroplasticity':         ['Grin2b'],
}

# Not detected in any fraction:
NOT_DETECTED = ['Slc16a3', 'Slc16a7', 'Acss1', 'Pdk4', 'Nfe2l2', 'Keap1']

# Ordered gene list for rows
GENE_ORDER = []
for g_list in GENE_GROUPS.values():
    GENE_ORDER.extend(g_list)

# ── Fraction / condition config ────────────────────────────────────────────────
FRACTIONS = {
    'Chrom': 'Chromatin',
    'SN':    'Soluble nuclear',
    'Cyto':  'Cytosol',
    'Memb':  'Membrane',
}

COND_COLS = {
    'Intox': ('Fold change',   'Corrected'),
    'AW':    ('Fold change.1', 'Corrected.1'),
    'PA':    ('Fold change.2', 'Corrected.2'),
}

FRAC_COLORS = {
    'Chrom': '#C01E42',
    'SN':    '#1A5FA0',
    'Cyto':  '#E8A020',
    'Memb':  '#3AAA50',
}

COND_LABELS = {
    'Intox': 'Intoxication',
    'AW':    'Acute\nWithdrawal',
    'PA':    'Protracted\nAbstinence',
}

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# ── Extract data ───────────────────────────────────────────────────────────────
print('Loading data...')
rows = []
for frac_key, sheet_name in FRACTIONS.items():
    df = pd.read_excel(DATA_FILE, sheet_name=sheet_name)
    df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()

    # Apply Keep + Review filter for Chromatin only
    if frac_key == 'Chrom' and 'Filter' in df.columns:
        df = df[df['Filter'].isin(['Keep', 'Review'])]
        print(f'  Chromatin after Keep+Review filter: {len(df)} proteins')

    for gene in GENE_ORDER:
        match = df[df['Gene symbol'].str.lower() == gene.lower()]
        if match.empty:
            continue
        row_data = match.iloc[0]
        for cond_key, (fc_col, corr_col) in COND_COLS.items():
            fc_val   = pd.to_numeric(row_data.get(fc_col,   np.nan), errors='coerce')
            corr_val = pd.to_numeric(row_data.get(corr_col, np.nan), errors='coerce')
            sig      = (not np.isnan(fc_val)) and (not np.isnan(corr_val)) and \
                       (abs(fc_val) > FC_THRESH) and (corr_val > CORR_THRESH)
            rows.append({
                'Gene':      gene,
                'Fraction':  frac_key,
                'Condition': cond_key,
                'FC':        fc_val,
                'Corrected': corr_val,
                'Sig':       sig,
            })

long_df = pd.DataFrame(rows)

# ── Save Excel table ───────────────────────────────────────────────────────────
# Wide format: genes as rows, fraction×condition as columns (FC + Corr)
records = {}
for _, r in long_df.iterrows():
    g = r['Gene']
    if g not in records:
        records[g] = {'Gene': g}
    col_fc   = f"{r['Fraction']}_{r['Condition']}_FC"
    col_corr = f"{r['Fraction']}_{r['Condition']}_Corr"
    col_sig  = f"{r['Fraction']}_{r['Condition']}_Sig"
    records[g][col_fc]   = round(r['FC'],   4) if not np.isnan(r['FC'])   else np.nan
    records[g][col_corr] = round(r['Corrected'], 4) if not np.isnan(r['Corrected']) else np.nan
    records[g][col_sig]  = r['Sig']

wide_tbl = pd.DataFrame(list(records.values()))
# Put gene order first
wide_tbl['_order'] = wide_tbl['Gene'].map({g: i for i, g in enumerate(GENE_ORDER)})
wide_tbl = wide_tbl.sort_values('_order').drop(columns='_order')

# Add a "not detected" note sheet
not_det_df = pd.DataFrame({'Gene': NOT_DETECTED,
                            'Note': ['Not detected in dataset'] * len(NOT_DETECTED)})

with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    wide_tbl.to_excel(writer, sheet_name='FC_and_Corrected', index=False)
    not_det_df.to_excel(writer, sheet_name='Not_Detected', index=False)
print(f'Saved: {OUT_XLSX}')

# ── Build heatmap matrices ─────────────────────────────────────────────────────
# Columns: Intox×4fracs | AW×4fracs | PA×4fracs
CONDITIONS = ['Intox', 'AW', 'PA']
FRAC_KEYS  = ['Chrom', 'SN', 'Cyto', 'Memb']
COL_ORDER  = [(c, f) for c in CONDITIONS for f in FRAC_KEYS]

n_genes = len(GENE_ORDER)
n_cols  = len(COL_ORDER)

fc_mat   = np.full((n_genes, n_cols), np.nan)
sig_mat  = np.zeros((n_genes, n_cols), dtype=bool)

for gi, gene in enumerate(GENE_ORDER):
    for ci, (cond, frac) in enumerate(COL_ORDER):
        sub = long_df[(long_df['Gene'] == gene) &
                      (long_df['Condition'] == cond) &
                      (long_df['Fraction']  == frac)]
        if not sub.empty:
            fc_mat[gi, ci]  = sub.iloc[0]['FC']
            sig_mat[gi, ci] = sub.iloc[0]['Sig']

# ── Figure ─────────────────────────────────────────────────────────────────────
print('Building figure...')

CELL_W  = 0.65   # inches per column
CELL_H  = 0.55   # inches per gene row
LEFT    = 2.2    # left margin for gene labels + group labels
RIGHT   = 2.0    # right margin for colorbar + legend
TOP     = 1.8    # top margin for column headers
BOT     = 0.5

FIG_W = LEFT + n_cols * CELL_W + RIGHT
FIG_H = TOP  + n_genes * CELL_H + BOT

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

# Main axes for heatmap
ax_l   = LEFT  / FIG_W
ax_b   = BOT   / FIG_H
ax_w   = (n_cols * CELL_W) / FIG_W
ax_h   = (n_genes * CELL_H) / FIG_H

ax = fig.add_axes([ax_l, ax_b, ax_w, ax_h])

# Colormap
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'VP', ['#1A5FA0', '#FFFFFF', '#C01E42'])
norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)

im = ax.imshow(fc_mat, aspect='auto', cmap=cmap, norm=norm,
               interpolation='nearest')

# ── Significance markers ───────────────────────────────────────────────────────
for gi in range(n_genes):
    for ci in range(n_cols):
        if sig_mat[gi, ci]:
            ax.text(ci, gi, '*', ha='center', va='center',
                    fontsize=14, color='#111111', fontweight='bold')

# ── Grid lines between condition groups ───────────────────────────────────────
group_bounds = [4, 8]   # after Intox (4 cols), after AW (8 cols)
for x in group_bounds:
    ax.axvline(x - 0.5, color='#333333', lw=2.0)

# Fine grid
for x in range(n_cols - 1):
    ax.axvline(x + 0.5, color='#AAAAAA', lw=0.4, zorder=2)
for y in range(n_genes - 1):
    ax.axhline(y + 0.5, color='#AAAAAA', lw=0.4, zorder=2)

# Group separator lines between gene groups
group_sep = []
count = 0
for g_list in GENE_GROUPS.values():
    count += len(g_list)
    if count < n_genes:
        group_sep.append(count - 0.5)
for y in group_sep:
    ax.axhline(y, color='#444444', lw=1.5, zorder=3)

# ── Axes ticks & labels ────────────────────────────────────────────────────────
ax.set_xticks(range(n_cols))
frac_labels_row = [FRAC_KEYS[ci % 4] for ci in range(n_cols)]
ax.set_xticklabels(frac_labels_row, fontsize=9, fontfamily='Arial', rotation=45, ha='right')

ax.set_yticks(range(n_genes))
ax.set_yticklabels(GENE_ORDER, fontsize=11, fontfamily='Arial', fontstyle='italic')

ax.tick_params(axis='both', which='both', length=0, pad=4)

# ── Condition headers above heatmap ───────────────────────────────────────────
header_y = 1.0 + (0.5 / (n_genes * CELL_H))
for ci_cond, cond_key in enumerate(CONDITIONS):
    x_center = (ci_cond * 4 + 1.5) / n_cols
    ax.text(x_center, 1.045,
            COND_LABELS[cond_key],
            transform=ax.transAxes,
            ha='center', va='bottom',
            fontsize=12, fontweight='bold',
            fontfamily='Arial', color='#111111')
    # Bracket line
    x0 = (ci_cond * 4 + 0.05) / n_cols
    x1 = (ci_cond * 4 + 3.95) / n_cols
    ax.annotate('', xy=(x0, 1.025), xytext=(x1, 1.025),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='-', color='#444444', lw=1.5))

# ── Fraction color bar at top (thin colored strip) ────────────────────────────
strip_h = 0.018
strip_b = 1.033
for ci, (cond, frac) in enumerate(COL_ORDER):
    rect = matplotlib.patches.FancyBboxPatch(
        ((ci - 0.48) / n_cols, strip_b - strip_h),
        (0.96 / n_cols), strip_h,
        boxstyle='square,pad=0',
        transform=ax.transAxes,
        facecolor=FRAC_COLORS[frac], edgecolor='none',
        clip_on=False, zorder=5)
    ax.add_patch(rect)

# ── Group labels (left of heatmap) ────────────────────────────────────────────
# get_yaxis_transform(): x in axes fraction, y in DATA coordinates (0..n_genes-1)
count = 0
for grp_label, g_list in GENE_GROUPS.items():
    mid = count + (len(g_list) - 1) / 2   # data-coordinate midpoint
    ax.text(-0.16, mid,
            grp_label,
            transform=ax.get_yaxis_transform(),
            ha='right', va='center',
            fontsize=8.5, fontfamily='Arial', color='#555555',
            multialignment='center')
    # Bracket spanning the group
    y0 = count - 0.45
    y1 = count + len(g_list) - 0.55
    ax.annotate('', xy=(-0.07, y0), xytext=(-0.07, y1),
                xycoords=ax.get_yaxis_transform(),
                arrowprops=dict(arrowstyle='-', color='#888888', lw=1.2))
    count += len(g_list)

# ── Colorbar ──────────────────────────────────────────────────────────────────
cbar_l = ax_l + ax_w + 0.06
cbar_b = ax_b + ax_h * 0.2
cbar_w = 0.018
cbar_h = ax_h * 0.55

cax = fig.add_axes([cbar_l, cbar_b, cbar_w, cbar_h])
cb  = fig.colorbar(im, cax=cax, extend='both')
cb.set_label('Log$_2$ FC vs Naive', fontsize=10, fontfamily='Arial', labelpad=8)
cb.ax.tick_params(labelsize=9)
cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])

# ── Legend (fraction colors + sig marker) ─────────────────────────────────────
leg_x = cbar_l + 0.09
leg_y = ax_b + ax_h * 0.95

handles = []
for fk in FRAC_KEYS:
    short_names = {'Chrom': 'Chromatin', 'SN': 'Soluble Nuclear',
                   'Cyto': 'Cytosol', 'Memb': 'Membrane'}
    handles.append(mpatches.Patch(facecolor=FRAC_COLORS[fk],
                                  label=short_names[fk],
                                  edgecolor='#333333', linewidth=0.5))

fig.legend(handles=handles,
           loc='upper right',
           bbox_to_anchor=(0.98, 0.92),
           fontsize=9, framealpha=0.9,
           edgecolor='#CCCCCC', title='Fraction',
           title_fontsize=9)

# Sig annotation note
fig.text(cbar_l + cbar_w / 2, ax_b - 0.04,
         '* |FC|>0.5\n  corr.p>3.3',
         ha='center', va='top', fontsize=8,
         fontfamily='Arial', color='#333333')

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(0.5, 0.98,
         'Candidate Genes — Log$_2$ FC vs Naive Across All Fractions & Conditions',
         ha='center', va='top', fontsize=14, fontweight='bold', fontfamily='Arial')

fig.text(0.5, 0.008,
         f'Detected: {len(GENE_ORDER)} genes  |  Not detected in dataset: '
         f'{", ".join(NOT_DETECTED)}  |  AW-M-3 excluded',
         ha='center', va='bottom', fontsize=8, color='#666666',
         style='italic', fontfamily='Arial')

# ── Save ──────────────────────────────────────────────────────────────────────
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
