"""
Metabolic gene heatmap — FC across all fractions & conditions
=============================================================
Rows: all genes found in ≥1 fraction.
Columns: 12 (4 fractions × 3 conditions), FC (log2).
Grey = not detected. Dot = significant (|FC|>0.5, corr>3.3).
Dynamic VMAX. Two outputs: AllFound and SigOnly.

Outputs: MetabolicGenes_Heatmap_AllFound.pdf/.png
         MetabolicGenes_Heatmap_SigOnly.pdf/.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')

FC_THRESH   = 0.5
CORR_THRESH = 3.3

GENE_LIST = [
    ('ACSS2',           'Acss2'),
    ('ACSS1',           'Acss1'),
    ('ACLY',            'Acly'),
    ('ACACA',           'Acaca'),
    ('SLC16A1 (MCT1)',  'Slc16a1'),
    ('SLC16A7 (MCT2)',  'Slc16a7'),
    ('SLC16A3 (MCT4)',  'Slc16a3'),
    ('CS',              'Cs'),
    ('ACO2',            'Aco2'),
    ('IDH2',            'Idh2'),
    ('DLAT',            'Dlat'),
    ('PDHA1',           'Pdha1'),
    ('SUCLG1',          'Suclg1'),
    ('SDHA',            'Sdha'),
    ('SLC2A1 (GLUT1)',  'Slc2a1'),
    ('SLC2A3 (GLUT3)',  'Slc2a3'),
    ('HK1',             'Hk1'),
    ('PFKM',            'Pfkm'),
    ('ALDOA',           'Aldoa'),
    ('GAPDH',           'Gapdh'),
    ('PKM',             'Pkm'),
    ('LDHA',            'Ldha'),
    ('LDHB',            'Ldhb'),
    ('ADH1',            'Adh1'),
    ('CYP2E1',          'Cyp2e1'),
    ('CAT',             'Cat'),
    ('ADH5',            'Adh5'),
    ('ALDH2',           'Aldh2'),
    ('ALDH1A1',         'Aldh1a1'),
    ('SIRT1',           'Sirt1'),
    ('SIRT3',           'Sirt3'),
    ('NAMPT',           'Nampt'),
    ('PARP1',           'Parp1'),
    ('NFE2L2 (NRF2)',   'Nfe2l2'),
    ('SOD2',            'Sod2'),
    ('GPX4',            'Gpx4'),
    ('PRDX3',           'Prdx3'),
    ('TXN',             'Txn'),
    ('GCLC',            'Gclc'),
    ('EP300 (p300)',     'Ep300'),
    ('CREBBP (CBP)',     'Crebbp'),
    ('KAT2A',           'Kat2a'),
    ('HDAC1',           'Hdac1'),
    ('HDAC2',           'Hdac2'),
    ('SIRT6',           'Sirt6'),
    ('GFAP',            'Gfap'),
    ('GLUL',            'Glul'),
    ('SLC1A2 (EAAT2)',  'Slc1a2'),
    ('SLC1A3 (EAAT1)',  'Slc1a3'),
    ('GOT1',            'Got1'),
    ('SHMT1',           'Shmt1'),
    ('DLD',             'Dld'),
    ('FH',              'Fh'),
    ('MDH2',            'Mdh2'),
    ('MDH1',            'Mdh1'),
    ('GOT2',            'Got2'),
    ('OGDH',            'Ogdh'),
    ('ATP5A1',          'Atp5a1'),
    ('NDUFS1',          'Ndufs1'),
    ('UQCRC1',          'Uqcrc1'),
    ('COX4I1',          'Cox4i1'),
    ('TFAM',            'Tfam'),
    ('PPARGC1A (PGC1a)','Ppargc1a'),
    ('HIF1A',           'Hif1a'),
    ('CREB1',           'Creb1'),
    ('CAMK2A',          'Camk2a'),
    ('CAMK2B',          'Camk2b'),
    ('PRKAA1 (AMPK)',   'Prkaa1'),
    ('PRKAA2 (AMPK)',   'Prkaa2'),
    ('MTOR',            'Mtor'),
    ('RPTOR',           'Rptor'),
    ('RICTOR',          'Rictor'),
    ('FOXO3',           'Foxo3'),
    ('HSPA9',           'Hspa9'),
    ('VDAC1',           'Vdac1'),
    ('TOMM20',          'Tomm20'),
    ('CPT1A',           'Cpt1a'),
    ('CPT2',            'Cpt2'),
    ('ACADM',           'Acadm'),
    ('HADHA',           'Hadha'),
    ('HADHB',           'Hadhb'),
    ('G6PD',            'G6pd'),
    ('PGD',             'Pgd'),
    ('ME1',             'Me1'),
    ('ME2',             'Me2'),
    ('GLUD1',           'Glud1'),
    ('GLS',             'Gls'),
    ('PC',              'Pc'),
    ('PDC',             'Pdc'),
    ('PDK1',            'Pdk1'),
    ('PDK2',            'Pdk2'),
    ('FASN',            'Fasn'),
    ('SREBF1',          'Srebf1'),
    ('HMOX1',           'Hmox1'),
    ('NOS2',            'Nos2'),
    ('ARG1',            'Arg1'),
    ('IL1B',            'Il1b'),
    ('TNF',             'Tnf'),
    ('RELA',            'Rela'),
    ('NFKB1',           'Nfkb1'),
    ('ATF4',            'Atf4'),
    ('DDIT3 (CHOP)',    'Ddit3'),
    ('XBP1',            'Xbp1'),
    ('EIF2AK3 (PERK)',  'Eif2ak3'),
    ('IRE1 (ERN1)',      'Ern1'),
    ('ATG5',            'Atg5'),
    ('ATG7',            'Atg7'),
    ('MAP1LC3B (LC3B)', 'Map1lc3b'),
    ('SQSTM1 (p62)',    'Sqstm1'),
    ('BNIP3',           'Bnip3'),
    ('PINK1',           'Pink1'),
    ('PRKN (Parkin)',   'Prkn'),
    ('PDHB',            'Pdhb'),
    ('PDHX',            'Pdhx'),
]

FRACTIONS = [
    ('Membrane',        'Membrane', 'Memb'),
    ('Cytosol',         'Cytosol',  'Cyto'),
    ('Chromatin',       'Chromatin','Chrom'),
    ('Soluble Nuclear', 'Soluble nuclear', 'SN'),
]

COND_COLS = {
    'Intox': ('Fold change',   'Corrected'),
    'AW':    ('Fold change.1', 'Corrected.1'),
    'PA':    ('Fold change.2', 'Corrected.2'),
}
CONDITIONS = ['Intox', 'AW', 'PA']
COND_COLORS = ['#7B9DC4', '#C4827B', '#8EB87D']   # blue/red/green — column header tints

# ── Load sheets ────────────────────────────────────────────────────────────────
sheet_cache = {}
def get_sheet(sheet_name):
    if sheet_name not in sheet_cache:
        df = pd.read_excel(DATA_FILE, sheet_name=sheet_name)
        df['Gene symbol'] = df['Gene symbol'].astype(str).str.strip()
        if 'Filter' in df.columns:
            df = df[df['Filter'].isin(['Keep', 'Review'])]
        sheet_cache[sheet_name] = df
    return sheet_cache[sheet_name]

# ── Collect data ───────────────────────────────────────────────────────────────
# fc_matrix[gene_i, col_i]  = FC value or NaN
# sig_matrix[gene_i, col_i] = True/False
# not_detected[gene_i, col_i] = True if gene absent from that fraction entirely

n_genes = len(GENE_LIST)
n_cols  = len(FRACTIONS) * len(CONDITIONS)   # 12

fc_matrix  = np.full((n_genes, n_cols), np.nan)
sig_matrix = np.zeros((n_genes, n_cols), dtype=bool)
missing    = np.zeros((n_genes, n_cols), dtype=bool)   # not in that fraction at all

gene_labels = [d[0] for d in GENE_LIST]

for gi, (display, rat_sym) in enumerate(GENE_LIST):
    for fi, (frac_label, sheet_name, _) in enumerate(FRACTIONS):
        df    = get_sheet(sheet_name)
        match = df[df['Gene symbol'].str.lower() == rat_sym.lower()]
        base_col = fi * len(CONDITIONS)

        if match.empty:
            for ci in range(len(CONDITIONS)):
                missing[gi, base_col + ci] = True
        else:
            r = match.iloc[0]
            for ci, cond in enumerate(CONDITIONS):
                fc_col, corr_col = COND_COLS[cond]
                fc   = pd.to_numeric(r.get(fc_col,   np.nan), errors='coerce')
                corr = pd.to_numeric(r.get(corr_col, np.nan), errors='coerce')
                col_i = base_col + ci
                fc_matrix[gi, col_i] = fc
                if (not pd.isna(fc) and not pd.isna(corr)
                        and corr > CORR_THRESH and abs(fc) > FC_THRESH):
                    sig_matrix[gi, col_i] = True

# ── Filter helpers ─────────────────────────────────────────────────────────────
found_mask = ~missing.all(axis=1)          # at least one fraction detected
sig_mask   = sig_matrix.any(axis=1)       # significant in ≥1 col

def draw_heatmap(gene_idx, out_pdf, out_png, title):
    labels  = [gene_labels[i] for i in gene_idx]
    fc_sub  = fc_matrix[gene_idx, :]
    sig_sub = sig_matrix[gene_idx, :]
    mis_sub = missing[gene_idx, :]

    n = len(gene_idx)

    # Dynamic VMAX
    valid = fc_sub[~np.isnan(fc_sub)]
    vmax  = max(float(np.ceil(np.abs(valid).max())), 1.0) if len(valid) else 3.0

    # Colormap: red-white-blue (up=red, down=blue)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color='#E0E0E0')   # grey for NaN (not detected)

    # Mask: NaN for not-detected cells so they show grey
    plot_data = fc_sub.copy()
    plot_data[mis_sub] = np.nan

    # Figure sizing
    row_h   = 0.32   # inches per gene row
    fig_h   = max(8, n * row_h + 3.0)
    fig_w   = 10.0

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    im = ax.imshow(plot_data, aspect='auto', cmap=cmap,
                   vmin=-vmax, vmax=vmax,
                   interpolation='nearest')

    # Significance dots
    for ri in range(n):
        for ci in range(n_cols):
            if sig_sub[ri, ci]:
                ax.plot(ci, ri, 'o', color='#111111',
                        markersize=3.5, markeredgewidth=0, zorder=3)

    # Grid lines between fractions (every 3 cols)
    for x in [2.5, 5.5, 8.5]:
        ax.axvline(x, color='#333333', linewidth=1.8, zorder=4)

    # Thin grid between all cells
    for x in np.arange(-0.5, n_cols, 1):
        ax.axvline(x, color='white', linewidth=0.4, zorder=2)
    for y in np.arange(-0.5, n, 1):
        ax.axhline(y, color='white', linewidth=0.4, zorder=2)

    # X-axis: condition labels, colored tick labels
    col_labels = []
    for frac_label, _, _ in FRACTIONS:
        for cond in CONDITIONS:
            col_labels.append(cond)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=8, fontfamily='Arial', rotation=45, ha='right')
    for tick, ci in zip(ax.get_xticklabels(), range(n_cols)):
        tick.set_color(COND_COLORS[ci % 3])
        tick.set_fontweight('bold')

    # Fraction labels above
    for fi, (frac_label, _, _) in enumerate(FRACTIONS):
        mid = fi * 3 + 1
        ax.text(mid, -1.8, frac_label, ha='center', va='bottom',
                fontsize=9, fontfamily='Arial', fontweight='bold', color='#222222',
                transform=ax.get_xaxis_transform())

    # Y-axis: gene labels
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=7.5, fontfamily='Arial', fontstyle='italic')

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.tick_params(axis='x', length=0, pad=2)
    ax.tick_params(axis='y', length=0, pad=4)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02, aspect=20)
    cbar.set_label('Log2 FC vs Naïve', fontsize=9, fontfamily='Arial')
    cbar.ax.tick_params(labelsize=8)
    tick_step = 1 if vmax <= 4 else 2
    cbar.set_ticks(range(int(-vmax), int(vmax) + 1, tick_step))

    title_clean = title.replace('₂', '2')
    ax.set_title(title_clean, fontsize=13, fontweight='bold',
                 fontfamily='Arial', pad=28)

    fig.text(0.5, 0.005,
             '● significant (|FC| > 0.5, corr. −log2 p > 3.3)  |  '
             'Grey = not detected  |  AW-M-3 excluded  |  Chromatin: Keep+Review',
             ha='center', va='bottom', fontsize=7,
             color='#888888', style='italic', fontfamily='Arial')

    plt.subplots_adjust(left=0.22, right=0.92, top=0.93, bottom=0.08)

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_pdf}')
    print(f'Saved: {out_png}')

# ── All found ──────────────────────────────────────────────────────────────────
all_idx = np.where(found_mask)[0]
draw_heatmap(
    all_idx,
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_Heatmap_AllFound.pdf'),
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_Heatmap_AllFound.png'),
    'Metabolic Candidate Genes — Log₂ FC Across Fractions & Conditions',
)

# ── Sig only ───────────────────────────────────────────────────────────────────
sig_idx = np.where(sig_mask)[0]
draw_heatmap(
    sig_idx,
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_Heatmap_SigOnly.pdf'),
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_Heatmap_SigOnly.png'),
    'Metabolic Candidate Genes — Significant Hits Only',
)

print(f'\nAll found: {len(all_idx)} genes')
print(f'Sig only:  {len(sig_idx)} genes')
