"""
Union Translocation Proteins — All-Fraction Table & Plot
=========================================================
Takes the 746 union translocation proteins and looks up their
FC (log2) and corrected p-value in ALL 4 fractions
(Chromatin, Soluble Nuclear, Cytosol, Membrane) for all 3
conditions (Intoxication, Acute Withdrawal, Protracted Abstinence).

Outputs:
  Translocation_Union_AllFractions_table.xlsx
  Translocation_Union_AllFractions_plot.pdf / .png
"""

import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..', '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
UNION_FILE = os.path.join(SCRIPT_DIR, 'Translocation_Union_stats.xlsx')
OUT_XLSX   = os.path.join(SCRIPT_DIR, 'Translocation_Union_AllFractions_table.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'Translocation_Union_AllFractions_plot.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'Translocation_Union_AllFractions_plot.png')

# ── Fraction config ────────────────────────────────────────────────────────────
# Each sheet's pre-computed columns:
#   Intox:  'Fold change',   'Corrected'
#   AW:     'Fold change.1', 'Corrected.1'
#   PA:     'Fold change.2', 'Corrected.2'
FRACTIONS = {
    'Chromatin':      'Chromatin',
    'SolubleNuclear': 'Soluble nuclear',
    'Cytosol':        'Cytosol',
    'Membrane':       'Membrane',
}

COND_COLS = {
    'Intox': ('Fold change',   'Corrected'),
    'AW':    ('Fold change.1', 'Corrected.1'),
    'PA':    ('Fold change.2', 'Corrected.2'),
}

FRAC_COLORS = {
    'Chromatin':      '#C01E42',   # dark red
    'SolubleNuclear': '#1A5FA0',   # dark blue
    'Cytosol':        '#E8A020',   # amber
    'Membrane':       '#3AAA50',   # green
}

FRAC_LABELS = {
    'Chromatin':      'Chromatin',
    'SolubleNuclear': 'Soluble Nuclear',
    'Cytosol':        'Cytosol',
    'Membrane':       'Membrane',
}

COND_LABELS = {
    'Intox': 'Intoxication',
    'AW':    'Acute Withdrawal',
    'PA':    'Protracted Abstinence',
}

# ── Load union proteins ────────────────────────────────────────────────────────
print('Loading union proteins...')
union_df = pd.read_excel(UNION_FILE)
union_746 = union_df[
    union_df['Sig_Intox'] | union_df['Sig_AW'] | union_df['Sig_PA']
].copy()
print(f'Union proteins: {len(union_746)}')
union_genes = set(union_746['Gene'].astype(str).str.strip())

# ── Load each fraction and extract FC + Corrected for union proteins ───────────
print('Loading fraction sheets...')
all_rows = []

for frac_key, sheet_name in FRACTIONS.items():
    df = pd.read_excel(DATA_FILE, sheet_name=sheet_name)
    df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
    df_frac = df[df['Gene symbol'].isin(union_genes)].copy()
    print(f'  {frac_key}: {len(df_frac)} of {len(union_746)} proteins found')

    for cond_key, (fc_col, corr_col) in COND_COLS.items():
        for _, row in df_frac.iterrows():
            fc_val   = pd.to_numeric(row.get(fc_col,   np.nan), errors='coerce')
            corr_val = pd.to_numeric(row.get(corr_col, np.nan), errors='coerce')
            all_rows.append({
                'Gene':      row['Gene symbol'],
                'Fraction':  frac_key,
                'Condition': cond_key,
                'FC':        fc_val,
                'Corrected': corr_val,
            })

long_df = pd.DataFrame(all_rows)

# ── Build wide table ───────────────────────────────────────────────────────────
print('Building wide table...')
pivot_fc   = long_df.pivot_table(index='Gene', columns=['Fraction', 'Condition'],
                                  values='FC',        aggfunc='first')
pivot_corr = long_df.pivot_table(index='Gene', columns=['Fraction', 'Condition'],
                                  values='Corrected', aggfunc='first')

# Flatten column names
pivot_fc.columns   = [f'{f}_{c}_FC'   for f, c in pivot_fc.columns]
pivot_corr.columns = [f'{f}_{c}_Corr' for f, c in pivot_corr.columns]

wide = pd.concat([pivot_fc, pivot_corr], axis=1).reset_index()

# Add translocation metadata from union
meta_cols = ['Gene', 'Direction_Intox', 'Sig_Intox',
                     'Direction_AW',    'Sig_AW',
                     'Direction_PA',    'Sig_PA']
wide = union_746[meta_cols].rename(columns={'Gene': 'Gene'}).merge(
    wide, on='Gene', how='right')

# Sort by AW Chromatin FC descending
aw_ch_col = 'Chromatin_AW_FC'
if aw_ch_col in wide.columns:
    wide = wide.sort_values(aw_ch_col, ascending=False, na_position='last')

wide.to_excel(OUT_XLSX, index=False)
print(f'Saved: {OUT_XLSX}')

# ── Plot ───────────────────────────────────────────────────────────────────────
print('Building plot...')

# For plot: use union proteins sorted by AW Chromatin FC
plot_genes = wide['Gene'].tolist()
n_genes    = len(plot_genes)
gene_idx   = {g: i for i, g in enumerate(plot_genes)}

# Figure size: very wide to fit all genes
FIG_W   = max(60, n_genes * 0.10)
FIG_H   = 18
MARKER  = 'o'
MS      = 3       # marker size
ALPHA   = 0.70
LW      = 0.0     # no edge

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
fig.suptitle(
    f'Translocation Union Proteins — FC vs Naïve Across All Fractions  (n = {n_genes})',
    fontsize=20, fontweight='bold', y=0.98, fontfamily='Arial')

gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.45,
                       top=0.93, bottom=0.10, left=0.03, right=0.985)

SHOW_XLABELS = n_genes <= 300  # only print gene names if ≤ 300 proteins

for ci, (cond_key, cond_label) in enumerate(COND_LABELS.items()):
    ax = fig.add_subplot(gs[ci])

    for frac_key in FRACTIONS:
        subset = long_df[
            (long_df['Condition'] == cond_key) &
            (long_df['Fraction']  == frac_key)
        ].copy()
        subset['x'] = subset['Gene'].map(gene_idx)
        subset = subset.dropna(subset=['x', 'FC'])

        ax.scatter(subset['x'].values, subset['FC'].values,
                   s=MS**2, color=FRAC_COLORS[frac_key],
                   alpha=ALPHA, linewidths=LW,
                   label=FRAC_LABELS[frac_key],
                   rasterized=True, zorder=3)

    # Reference lines
    ax.axhline(0,   color='#888888', lw=0.8, ls='-',  zorder=1)
    ax.axhline( 0.5, color='#CCCCCC', lw=0.6, ls='--', zorder=1)
    ax.axhline(-0.5, color='#CCCCCC', lw=0.6, ls='--', zorder=1)

    ax.set_xlim(-1, n_genes)
    ymax = max(3.5, long_df[long_df['Condition'] == cond_key]['FC'].abs().quantile(0.995) + 0.3)
    ax.set_ylim(-ymax, ymax)

    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    ax.set_ylabel('Log$_2$ FC vs Naive', fontsize=13, fontfamily='Arial')
    ax.set_title(cond_label, fontsize=15, fontweight='bold',
                 fontfamily='Arial', pad=6)
    ax.tick_params(labelsize=10)

    if SHOW_XLABELS:
        ax.set_xticks(range(n_genes))
        ax.set_xticklabels(plot_genes, rotation=90, fontsize=5,
                           fontfamily='Arial', ha='center')
    else:
        # Too many genes — show tick marks but no labels; annotate count
        ax.set_xticks([])
        ax.text(0.02, 0.96,
                f'{n_genes} proteins sorted by AW Chromatin FC (high → low)',
                transform=ax.transAxes, ha='left', va='top',
                fontsize=10, color='#555555', style='italic', fontfamily='Arial')

    # Legend on first panel only
    if ci == 0:
        handles = [mpatches.Patch(color=FRAC_COLORS[f], label=FRAC_LABELS[f])
                   for f in FRACTIONS]
        ax.legend(handles=handles, fontsize=11, loc='upper right',
                  framealpha=0.9, edgecolor='#CCCCCC', frameon=True,
                  borderpad=0.7, ncol=4)

# X-axis label on bottom panel only
fig.text(0.5, 0.02,
         'Gene (sorted by AW Chromatin FC, high → low)  |  '
         'AW-M-3 excluded  |  Chromatin: Keep+Review only',
         ha='center', va='bottom', fontsize=9, color='#666666',
         style='italic', fontfamily='Arial')

print(f'Saving PDF ({FIG_W:.0f} × {FIG_H} inches)...')
with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=150, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
print('Done.')
