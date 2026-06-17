"""
Simple candidate gene table — one compartment per gene (biological default)
Columns: Gene | Fraction | Intox FC | Intox Corr | AW FC | AW Corr | PA FC | PA Corr
Significant cells (|FC|>0.5 AND Corr>3.3) are highlighted.
"""

import os, pandas as pd, numpy as np
import matplotlib, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')
OUT_PDF    = os.path.join(SCRIPT_DIR, 'CandidateGenes_SimpleTable.pdf')
OUT_PNG    = os.path.join(SCRIPT_DIR, 'CandidateGenes_SimpleTable.png')

# Gene → home compartment(s) as a list of (sheet name, display label)
# Genes with multiple entries appear as separate rows
GENE_FRAC = [
    ('Slc16a1', 'Membrane',       'Membrane'),
    ('Slc2a1',  'Membrane',       'Membrane'),
    ('Slc2a3',  'Membrane',       'Membrane'),
    ('Pfkl',    'Cytosol',        'Cytosol'),
    ('Pfkm',    'Cytosol',        'Cytosol'),
    ('Pfkp',    'Cytosol',        'Cytosol'),
    ('Pdha1',   'Cytosol',        'Cytosol'),
    ('Pdk1',    'Cytosol',        'Cytosol'),
    ('Pdk2',    'Cytosol',        'Cytosol'),
    ('Pdk3',    'Cytosol',        'Cytosol'),
    ('Acss2',   'Cytosol',        'Cytosol'),
    ('Acss2',   'Chromatin',      'Chromatin'),
    ('Acly',    'Cytosol',        'Cytosol'),
    ('Adh5',    'Cytosol',        'Cytosol'),
    ('Aldh2',   'Cytosol',        'Cytosol'),
    ('Aldh2',   'Chromatin',      'Chromatin'),
    ('Cyp2e1',  'Membrane',       'Membrane'),
    ('Cat',     'Cytosol',        'Cytosol'),
    ('Grin2b',  'Membrane',       'Membrane'),
]

COND_COLS = {
    'Intox': ('Fold change',   'Corrected'),
    'AW':    ('Fold change.1', 'Corrected.1'),
    'PA':    ('Fold change.2', 'Corrected.2'),
}

CORR_THRESH = 3.3
FC_THRESH   = 0.5

# ── Pull data ──────────────────────────────────────────────────────────────────
sheet_cache = {}
def get_sheet(sheet_name):
    if sheet_name not in sheet_cache:
        df = pd.read_excel(DATA_FILE, sheet_name=sheet_name)
        df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
        if 'Filter' in df.columns:          # Chromatin only — keep+review
            df = df[df['Filter'].isin(['Keep', 'Review'])]
        sheet_cache[sheet_name] = df
    return sheet_cache[sheet_name]

rows = []
for gene, sheet, frac_label in GENE_FRAC:
    df = get_sheet(sheet)
    match = df[df['Gene symbol'].str.lower() == gene.lower()]
    entry = {'Gene': gene, 'Fraction': frac_label}
    for cond, (fc_col, corr_col) in COND_COLS.items():
        if match.empty:
            entry[f'{cond}_FC']   = np.nan
            entry[f'{cond}_Corr'] = np.nan
        else:
            r = match.iloc[0]
            entry[f'{cond}_FC']   = pd.to_numeric(r.get(fc_col,   np.nan), errors='coerce')
            entry[f'{cond}_Corr'] = pd.to_numeric(r.get(corr_col, np.nan), errors='coerce')
    rows.append(entry)

data = pd.DataFrame(rows)

# ── Draw table ─────────────────────────────────────────────────────────────────
CONDITIONS = ['Intox', 'AW', 'PA']
COND_NICE  = {'Intox': 'Intoxication', 'AW': 'Acute Withdrawal', 'PA': 'Protracted Abstinence'}

n_genes = len(data)

# Column layout: Gene | Fraction | [FC  Corr] × 3
COL_WIDTHS  = [1.2, 1.1, 0.7, 0.85, 0.7, 0.85, 0.7, 0.85]   # inches
ROW_H       = 0.32   # inches per row
HEADER_H    = 0.70   # inches for two-row header

FIG_W = sum(COL_WIDTHS) + 0.3
FIG_H = HEADER_H + n_genes * ROW_H + 0.5

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_xlim(0, sum(COL_WIDTHS))
ax.set_ylim(0, FIG_H - 0.3)
ax.axis('off')

# cumulative x positions
xs = [0]
for w in COL_WIDTHS:
    xs.append(xs[-1] + w)
col_centers = [(xs[i] + xs[i+1]) / 2 for i in range(len(COL_WIDTHS))]

TOTAL_H = FIG_H - 0.3
header_y  = TOTAL_H - HEADER_H
data_top  = header_y

# Colors
C_UP   = '#FADADD'   # light red
C_DOWN = '#DAE8FA'   # light blue
C_HEAD = '#2C2C2C'
C_ALT  = '#F7F7F7'   # alternating row
C_LINE = '#BBBBBB'

def sig(fc, corr):
    if np.isnan(fc) or np.isnan(corr): return None
    if fc >  FC_THRESH and corr > CORR_THRESH: return 'up'
    if fc < -FC_THRESH and corr > CORR_THRESH: return 'down'
    return None

def fmt(val):
    if np.isnan(val): return '—'
    return f'{val:.2f}'

# ── Header ─────────────────────────────────────────────────────────────────────
# Top row: condition span labels
cond_spans = [(2, 4, 'Intoxication'), (4, 6, 'Acute Withdrawal'), (6, 8, 'Protracted Abstinence')]
for c0, c1, label in cond_spans:
    cx = (xs[c0] + xs[c1]) / 2
    cy = TOTAL_H - HEADER_H * 0.28
    ax.text(cx, cy, label, ha='center', va='center',
            fontsize=9, fontweight='bold', fontfamily='Arial', color='white')
    rect = plt.Rectangle((xs[c0] + 0.02, TOTAL_H - HEADER_H * 0.55),
                          xs[c1] - xs[c0] - 0.04, HEADER_H * 0.45,
                          facecolor='#3A3A3A', edgecolor='none', zorder=0)
    ax.add_patch(rect)

# Second row: column labels
col_labels = ['Gene', 'Fraction', 'FC', 'Corr. p', 'FC', 'Corr. p', 'FC', 'Corr. p']
for ci, (cx, lbl) in enumerate(zip(col_centers, col_labels)):
    ax.text(cx, TOTAL_H - HEADER_H * 0.78, lbl,
            ha='center', va='center',
            fontsize=8.5, fontweight='bold', fontfamily='Arial', color='#111111')

# Header bottom border
ax.axhline(header_y, color='#333333', lw=1.5, xmin=0, xmax=1)
ax.axhline(TOTAL_H,  color='#333333', lw=1.5, xmin=0, xmax=1)

# ── Data rows ─────────────────────────────────────────────────────────────────
for ri, row in data.iterrows():
    row_y   = data_top - (ri + 1) * ROW_H
    row_mid = row_y + ROW_H / 2

    # Alternating background
    if ri % 2 == 1:
        bg = plt.Rectangle((0, row_y), sum(COL_WIDTHS), ROW_H,
                            facecolor=C_ALT, edgecolor='none', zorder=0)
        ax.add_patch(bg)

    # Gene (italic) and Fraction
    ax.text(col_centers[0], row_mid, row['Gene'],
            ha='center', va='center', fontsize=9,
            fontstyle='italic', fontfamily='Arial', color='#111111')
    ax.text(col_centers[1], row_mid, row['Fraction'],
            ha='center', va='center', fontsize=8.5,
            fontfamily='Arial', color='#444444')

    # FC + Corr per condition
    for ci, cond in enumerate(CONDITIONS):
        fc_val   = row[f'{cond}_FC']
        corr_val = row[f'{cond}_Corr']
        direction = sig(fc_val, corr_val)

        col_fc   = 2 + ci * 2
        col_corr = 3 + ci * 2

        # Cell highlight for significant values
        for col_i in [col_fc, col_corr]:
            if direction:
                color = C_UP if direction == 'up' else C_DOWN
                rect = plt.Rectangle((xs[col_i] + 0.02, row_y + 0.02),
                                      COL_WIDTHS[col_i] - 0.04, ROW_H - 0.04,
                                      facecolor=color, edgecolor='none', zorder=1)
                ax.add_patch(rect)

        fc_text   = fmt(fc_val)
        corr_text = fmt(corr_val)
        fw = 'bold' if direction else 'normal'

        ax.text(col_centers[col_fc], row_mid, fc_text,
                ha='center', va='center', fontsize=8.5,
                fontfamily='Arial', fontweight=fw, color='#111111', zorder=2)
        ax.text(col_centers[col_corr], row_mid, corr_text,
                ha='center', va='center', fontsize=8.5,
                fontfamily='Arial', fontweight=fw, color='#111111', zorder=2)

    # Row bottom line
    ax.axhline(row_y, color=C_LINE, lw=0.5)

# Bottom border
ax.axhline(data_top - n_genes * ROW_H, color='#333333', lw=1.5)

# Vertical dividers
for x in [xs[1], xs[2], xs[4], xs[6]]:
    ax.axvline(x, color='#888888', lw=0.8,
               ymin=(data_top - n_genes * ROW_H) / (FIG_H - 0.3),
               ymax=1.0)

# ── Legend ────────────────────────────────────────────────────────────────────
leg_y = data_top - n_genes * ROW_H - 0.22
ax.add_patch(plt.Rectangle((0.1, leg_y - 0.10), 0.20, 0.14,
             facecolor=C_UP, edgecolor='#999999', lw=0.5))
ax.text(0.35, leg_y - 0.03, 'Increased (|FC|>0.5, corr. p>3.3)',
        fontsize=7.5, va='center', fontfamily='Arial', color='#333333')

ax.add_patch(plt.Rectangle((3.2, leg_y - 0.10), 0.20, 0.14,
             facecolor=C_DOWN, edgecolor='#999999', lw=0.5))
ax.text(3.45, leg_y - 0.03, 'Decreased',
        fontsize=7.5, va='center', fontfamily='Arial', color='#333333')

ax.text(sum(COL_WIDTHS), leg_y - 0.03,
        'AW-M-3 excluded  |  Chromatin: Keep+Review only',
        fontsize=7, ha='right', va='center',
        fontfamily='Arial', color='#888888', style='italic')

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(0.5, 0.99,
         'Candidate Genes — FC vs Naive in Primary Compartment',
         ha='center', va='top', fontsize=12, fontweight='bold', fontfamily='Arial')

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
