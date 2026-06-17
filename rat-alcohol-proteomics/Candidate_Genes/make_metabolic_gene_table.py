"""
Metabolic candidate gene table — all 4 compartments per gene
=============================================================
Columns: Gene | [Membrane: Intox/AW/PA FC] | [Cytosol: ...] | [Chromatin: ...] | [SN: ...]
FC cells highlighted if |FC|>0.5 AND corr>3.3.
If gene not found in a compartment: grey dash.

Outputs (two versions):
  MetabolicGenes_AllFound.pdf/.png   — all genes found in ≥1 compartment
  MetabolicGenes_SigOnly.pdf/.png    — genes significant in ≥1 condition/compartment
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE  = os.path.join(SCRIPT_DIR, '..',
             'EDIT  Excluding AWM3 Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx')

FC_THRESH   = 0.5
CORR_THRESH = 3.3

# Human display name → rat gene symbol (first-letter cap, rest lower)
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
]

FRACTIONS = [
    ('Membrane',       'Membrane', 'Memb'),
    ('Cytosol',        'Cytosol',  'Cyto'),
    ('Chromatin',      'Chromatin','Chrom'),
    ('Soluble nuclear','SN',       'Nuc'),
]

COND_COLS = {
    'Intox': ('Fold change',   'Corrected'),
    'AW':    ('Fold change.1', 'Corrected.1'),
    'PA':    ('Fold change.2', 'Corrected.2'),
}
CONDITIONS = ['Intox', 'AW', 'PA']

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

def sig(fc, corr):
    if pd.isna(fc) or pd.isna(corr): return None
    if fc >  FC_THRESH and corr > CORR_THRESH: return 'up'
    if fc < -FC_THRESH and corr > CORR_THRESH: return 'down'
    return None

def fmt(val):
    if pd.isna(val): return '—'
    return f'{val:+.2f}'

# ── Build data ─────────────────────────────────────────────────────────────────
print('Searching all compartments...')
rows = []
for display, rat_sym in GENE_LIST:
    entry = {'Gene': display, 'RatSym': rat_sym}
    found_any = False
    sig_any   = False

    for sheet_name, frac_label, _ in FRACTIONS:
        df = get_sheet(sheet_name)
        match = df[df['Gene symbol'].str.lower() == rat_sym.lower()]
        in_frac = not match.empty

        for cond in CONDITIONS:
            fc_col, corr_col = COND_COLS[cond]
            if in_frac:
                r    = match.iloc[0]
                fc   = pd.to_numeric(r.get(fc_col,   np.nan), errors='coerce')
                corr = pd.to_numeric(r.get(corr_col, np.nan), errors='coerce')
                found_any = True
                if sig(fc, corr): sig_any = True
            else:
                fc, corr = np.nan, np.nan
            entry[f'{frac_label}_{cond}_FC']   = fc
            entry[f'{frac_label}_{cond}_Corr'] = corr
            entry[f'{frac_label}_{cond}_InFrac'] = in_frac

    entry['found_any'] = found_any
    entry['sig_any']   = sig_any
    rows.append(entry)

all_data = pd.DataFrame(rows)
found_data = all_data[all_data['found_any']].reset_index(drop=True)
sig_data   = all_data[all_data['sig_any']].reset_index(drop=True)

print(f'Total genes in list:   {len(all_data)}')
print(f'Found in ≥1 fraction:  {len(found_data)}')
print(f'Sig in ≥1 condition:   {len(sig_data)}')

# ── Draw table ─────────────────────────────────────────────────────────────────
C_UP   = '#FADADD'
C_DOWN = '#DAE8FA'
C_ALT  = '#F7F7F7'
C_LINE = '#BBBBBB'
C_MISS = '#F0F0F0'   # gene not in this fraction

FRAC_LABELS  = [f[1] for f in FRACTIONS]   # Membrane, Cytosol, Chromatin, SN
FRAC_COLORS  = ['#E8F4E8', '#FFF8E8', '#FFF0F0', '#E8F0FF']

# Column layout: Gene (1.1") + 4 fractions × 3 conditions × 0.55" each
GENE_W  = 1.4
COND_W  = 0.55
N_FRAC  = 4
N_COND  = 3
ROW_H   = 0.26
HEAD_H  = 0.72

def draw_table(data, out_pdf, out_png, title):
    n = len(data)
    if n == 0:
        print(f'No rows for {title}, skipping.')
        return

    col_widths = [GENE_W] + [COND_W] * (N_FRAC * N_COND)
    xs = [0]
    for w in col_widths: xs.append(xs[-1] + w)
    total_w = xs[-1]
    col_cx  = [(xs[i]+xs[i+1])/2 for i in range(len(col_widths))]

    FIG_W = total_w + 0.3
    FIG_H = HEAD_H + n * ROW_H + 0.45
    TOTAL_H = FIG_H - 0.3

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, TOTAL_H)
    ax.axis('off')

    header_y = TOTAL_H - HEAD_H
    data_top  = header_y

    # ── Header row 1: fraction labels ─────────────────────────────────────────
    for fi, (flab, fcol) in enumerate(zip(FRAC_LABELS, FRAC_COLORS)):
        c0 = xs[1 + fi * N_COND]
        c1 = xs[1 + (fi+1) * N_COND]
        cx = (c0 + c1) / 2
        rect = plt.Rectangle((c0 + 0.02, TOTAL_H - HEAD_H * 0.52),
                              c1 - c0 - 0.04, HEAD_H * 0.42,
                              facecolor='#3A3A3A', edgecolor='none', zorder=0)
        ax.add_patch(rect)
        ax.text(cx, TOTAL_H - HEAD_H * 0.31, flab,
                ha='center', va='center', fontsize=8.5,
                fontweight='bold', fontfamily='Arial', color='white')

    # ── Header row 2: condition labels ────────────────────────────────────────
    cond_short = {'Intox': 'Intox', 'AW': 'AW', 'PA': 'PA'}
    col_idx = 1
    for fi in range(N_FRAC):
        for cond in CONDITIONS:
            cx = col_cx[col_idx]
            ax.text(cx, TOTAL_H - HEAD_H * 0.79, cond_short[cond],
                    ha='center', va='center', fontsize=7.5,
                    fontweight='bold', fontfamily='Arial', color='#222222')
            col_idx += 1

    # Gene header
    ax.text(col_cx[0], TOTAL_H - HEAD_H * 0.55, 'Gene',
            ha='center', va='center', fontsize=8.5,
            fontweight='bold', fontfamily='Arial', color='#111111')

    ax.axhline(TOTAL_H,  color='#333333', lw=1.5)
    ax.axhline(header_y, color='#333333', lw=1.5)

    # ── Data rows ─────────────────────────────────────────────────────────────
    for ri, row in data.iterrows():
        row_y   = data_top - (ri + 1) * ROW_H
        row_mid = row_y + ROW_H / 2

        if ri % 2 == 1:
            ax.add_patch(plt.Rectangle((0, row_y), total_w, ROW_H,
                         facecolor=C_ALT, edgecolor='none', zorder=0))

        ax.text(col_cx[0], row_mid, row['Gene'],
                ha='center', va='center', fontsize=7.5,
                fontstyle='italic', fontfamily='Arial', color='#111111')

        col_idx = 1
        for fi, (sheet_name, frac_label, _) in enumerate(FRACTIONS):
            for cond in CONDITIONS:
                fc      = row[f'{frac_label}_{cond}_FC']
                corr    = row[f'{frac_label}_{cond}_Corr']
                in_frac = row[f'{frac_label}_{cond}_InFrac']
                cx      = col_cx[col_idx]

                if not in_frac:
                    # Gene not in this fraction
                    ax.add_patch(plt.Rectangle(
                        (xs[col_idx]+0.01, row_y+0.01),
                        col_widths[col_idx]-0.02, ROW_H-0.02,
                        facecolor=C_MISS, edgecolor='none', zorder=1))
                    ax.text(cx, row_mid, '—', ha='center', va='center',
                            fontsize=7, fontfamily='Arial', color='#AAAAAA', zorder=2)
                else:
                    direction = sig(fc, corr)
                    if direction:
                        color = C_UP if direction == 'up' else C_DOWN
                        ax.add_patch(plt.Rectangle(
                            (xs[col_idx]+0.01, row_y+0.01),
                            col_widths[col_idx]-0.02, ROW_H-0.02,
                            facecolor=color, edgecolor='none', zorder=1))
                    fc_str = fmt(fc)
                    fw = 'bold' if direction else 'normal'
                    ax.text(cx, row_mid, fc_str, ha='center', va='center',
                            fontsize=7.5, fontfamily='Arial',
                            fontweight=fw, color='#111111', zorder=2)

                col_idx += 1

        ax.axhline(row_y, color=C_LINE, lw=0.4)

    ax.axhline(data_top - n * ROW_H, color='#333333', lw=1.5)

    # Vertical dividers between fractions
    for fi in range(N_FRAC + 1):
        x = xs[1 + fi * N_COND]
        ax.axvline(x, color='#888888', lw=0.8,
                   ymin=(data_top - n * ROW_H) / TOTAL_H,
                   ymax=header_y / TOTAL_H)

    # ── Legend ────────────────────────────────────────────────────────────────
    leg_y = data_top - n * ROW_H - 0.20
    ax.add_patch(plt.Rectangle((0.1, leg_y-0.10), 0.18, 0.13,
                 facecolor=C_UP, edgecolor='#999999', lw=0.5))
    ax.text(0.33, leg_y-0.035, 'Increased (|FC|>0.5, corr. p>3.3)',
            fontsize=7, va='center', fontfamily='Arial', color='#333333')
    ax.add_patch(plt.Rectangle((3.5, leg_y-0.10), 0.18, 0.13,
                 facecolor=C_DOWN, edgecolor='#999999', lw=0.5))
    ax.text(3.73, leg_y-0.035, 'Decreased',
            fontsize=7, va='center', fontfamily='Arial', color='#333333')
    ax.add_patch(plt.Rectangle((5.0, leg_y-0.10), 0.18, 0.13,
                 facecolor=C_MISS, edgecolor='#999999', lw=0.5))
    ax.text(5.23, leg_y-0.035, 'Not detected in fraction',
            fontsize=7, va='center', fontfamily='Arial', color='#333333')

    ax.text(total_w, leg_y-0.035,
            'FC = log2 fold change vs. naïve  |  AW-M-3 excluded  |  Chromatin: Keep+Review only',
            fontsize=6.5, ha='right', va='center',
            fontfamily='Arial', color='#888888', style='italic')

    # ── Title ─────────────────────────────────────────────────────────────────
    fig.text(0.5, 0.995, title,
             ha='center', va='top', fontsize=11,
             fontweight='bold', fontfamily='Arial')

    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {out_pdf}')
    print(f'Saved: {out_png}')

# ── Generate both versions ─────────────────────────────────────────────────────
draw_table(
    found_data,
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_AllFound.pdf'),
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_AllFound.png'),
    'Metabolic & Stress Candidate Genes — All 4 Compartments (Genes Detected in ≥1 Fraction)'
)

draw_table(
    sig_data,
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_SigOnly.pdf'),
    os.path.join(SCRIPT_DIR, 'MetabolicGenes_SigOnly.png'),
    'Metabolic & Stress Candidate Genes — Significant in ≥1 Condition (|FC|>0.5, corr. p>3.3)'
)

print('Done.')
