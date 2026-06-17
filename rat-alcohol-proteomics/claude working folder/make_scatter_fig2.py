import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.backends.backend_pdf import PdfPages

FILE    = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_PDF = 'Figure2_Scatter_Chrom_vs_Nuc_AW.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

C_CHROM = '#E8305A'   # pink  – shifted INTO Chromatin (above diagonal)
C_NUC   = '#2B7FD4'   # blue  – shifted INTO Nuclear   (below diagonal)
C_CO    = '#AAAAAA'   # gray  – co-regulated / not significant

COND_SUFFIXES = {
    'Naive':            ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Acute Withdrawal': ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
}

def calc_delta_and_stats(df, naive_cols, aw_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[aw_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        if len(n) >= 2 and len(c) >= 2:
            _, p = stats.ttest_ind(n, c, equal_var=False)
            p_vals.append(p)
        else:
            p_vals.append(np.nan)
    p_vals  = pd.Series(p_vals, index=df.index)
    valid   = p_vals.notna()
    ranks   = p_vals[valid].rank(ascending=False)
    corr    = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corr

def load_compartment(excel_name, prefix):
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    df  = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    naive_cols = [c for c in df.columns
                  if any(prefix + s in str(c) for s in COND_SUFFIXES['Naive'])]
    aw_cols    = [c for c in df.columns
                  if any(prefix + s in str(c) for s in COND_SUFFIXES['Acute Withdrawal'])]
    return df, naive_cols, aw_cols

# ── Load ──
nuc_df,   nuc_naive,   nuc_aw   = load_compartment('Soluble nuclear', 'Nuc_')
chrom_df, chrom_naive, chrom_aw = load_compartment('Chromatin',       'Chrom_')

nuc_fc,   nuc_corr   = calc_delta_and_stats(nuc_df,   nuc_naive,   nuc_aw)
chrom_fc, chrom_corr = calc_delta_and_stats(chrom_df, chrom_naive, chrom_aw)

nuc_sig_mask    = (nuc_corr   > CORR_THRESH) & (nuc_fc.abs()   > FC_THRESH)
chrom_sig_mask  = (chrom_corr > CORR_THRESH) & (chrom_fc.abs() > FC_THRESH)
nuc_sig_genes   = set(nuc_df.loc[nuc_sig_mask,    'Gene symbol'].astype(str))
chrom_sig_genes = set(chrom_df.loc[chrom_sig_mask, 'Gene symbol'].astype(str))
union_genes     = sorted(nuc_sig_genes | chrom_sig_genes)

nuc_delta_lk   = dict(zip(nuc_df['Gene symbol'].astype(str),   nuc_fc))
chrom_delta_lk = dict(zip(chrom_df['Gene symbol'].astype(str), chrom_fc))

nuc_lk   = (nuc_df.drop_duplicates('Gene symbol')
            .set_index(nuc_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))
chrom_lk = (chrom_df.drop_duplicates('Gene symbol')
            .set_index(chrom_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))

AW_NUC    = [c for c in nuc_df.columns
             if any('Nuc_'   + s in str(c) for s in COND_SUFFIXES['Acute Withdrawal'])]
AW_CHROM  = [c for c in chrom_df.columns
             if any('Chrom_' + s in str(c) for s in COND_SUFFIXES['Acute Withdrawal'])]
NAI_NUC   = [c for c in nuc_df.columns
             if any('Nuc_'   + s in str(c) for s in COND_SUFFIXES['Naive'])]
NAI_CHROM = [c for c in chrom_df.columns
             if any('Chrom_' + s in str(c) for s in COND_SUFFIXES['Naive'])]

def get_vals(lk, gene, cols):
    if gene not in lk.index: return np.full(len(cols), np.nan)
    row = lk.loc[gene]
    return np.array([float(pd.to_numeric(row[c], errors='coerce'))
                     if c in row.index else np.nan for c in cols])

rows = []
for gene in union_genes:
    dn = nuc_delta_lk.get(gene, np.nan)
    dc = chrom_delta_lk.get(gene, np.nan)
    if not (np.isfinite(dn) and np.isfinite(dc)):
        continue
    interaction = dc - dn          # > 0 → more change in Chrom
    chr_aw_v = get_vals(chrom_lk, gene, AW_CHROM)
    nuc_aw_v = get_vals(nuc_lk,   gene, AW_NUC)
    chr_n_v  = get_vals(chrom_lk, gene, NAI_CHROM)
    nuc_n_v  = get_vals(nuc_lk,   gene, NAI_NUC)
    cd = chr_aw_v[np.isfinite(chr_aw_v)] - np.nanmean(chr_n_v)
    nd = nuc_aw_v[np.isfinite(nuc_aw_v)] - np.nanmean(nuc_n_v)
    if len(cd) >= 2 and len(nd) >= 2:
        _, p = stats.ttest_ind(cd, nd, equal_var=False)
    else:
        p = np.nan
    rows.append({'Gene': gene, 'dNuc': dn, 'dChrom': dc,
                 'interaction': interaction, 'praw': p})

df = pd.DataFrame(rows)

# BH correction
valid = df['praw'].notna()
pv    = df.loc[valid, 'praw'].values
order = np.argsort(pv)
ranks = np.empty_like(order); ranks[order] = np.arange(1, len(pv)+1)
bh    = np.minimum(1.0, pv * len(pv) / ranks)
for i in range(len(bh)-2, -1, -1):
    bh[i] = min(bh[i], bh[i+1])
df.loc[valid, 'pBH'] = bh

# Classification: above diagonal = INTO Chromatin, below = INTO Nuclear
sig_into_chrom = (df['interaction'] > 0) & (df['pBH'] < 0.05)  # above y=x, significant
sig_into_nuc   = (df['interaction'] < 0) & (df['pBH'] < 0.05)  # below y=x, significant
ns_mask        = ~sig_into_chrom & ~sig_into_nuc

n_chrom = sig_into_chrom.sum()
n_nuc   = sig_into_nuc.sum()

# Pearson / Spearman on all proteins
r_p, p_p = pearsonr(df['dNuc'], df['dChrom'])
r_s, p_s = spearmanr(df['dNuc'], df['dChrom'])
print(f'INTO Chromatin (BH<0.05, above diagonal): {n_chrom}')
print(f'INTO Nuclear   (BH<0.05, below diagonal): {n_nuc}')
print(f'Pearson  r={r_p:.3f}  p={p_p:.2e}')
print(f'Spearman ρ={r_s:.3f}  p={p_s:.2e}')

# ── Plot ──
# Use actual data range + padding so no dots get clipped
all_x = df['dNuc'].dropna().values
all_y = df['dChrom'].dropna().values
pad   = 0.4
XLIM  = (min(all_x) - pad, max(all_x) + pad)
YLIM  = (min(all_y) - pad, max(all_y) + pad)
LIM   = max(abs(XLIM[0]), abs(XLIM[1]), abs(YLIM[0]), abs(YLIM[1]))  # for diagonal shading

with PdfPages(OUT_PDF) as pdf:
    fig, ax = plt.subplots(figsize=(8, 7.5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)

    # Diagonal shading: above y=x = pink (INTO Chrom), below = blue (INTO Nuc)
    tri_chrom = plt.Polygon([(-LIM, -LIM), (LIM, LIM), (-LIM, LIM)],
                             closed=True, color=C_CHROM, alpha=0.07, zorder=0)
    tri_nuc   = plt.Polygon([(-LIM, -LIM), (LIM, LIM), (LIM, -LIM)],
                             closed=True, color=C_NUC,   alpha=0.07, zorder=0)
    ax.add_patch(tri_chrom)
    ax.add_patch(tri_nuc)

    # Reference lines
    ax.axhline(0, color='#CCCCCC', linewidth=0.8, zorder=1)
    ax.axvline(0, color='#CCCCCC', linewidth=0.8, zorder=1)
    ax.plot([-LIM, LIM], [-LIM, LIM], color='#888888',
            linestyle='--', linewidth=1.4, zorder=2, label='y = x')

    # Scatter — gray behind, colored on top
    ax.scatter(df.loc[ns_mask,        'dNuc'], df.loc[ns_mask,        'dChrom'],
               s=12, color=C_CO,    alpha=0.30, linewidths=0, rasterized=True, zorder=3)
    ax.scatter(df.loc[sig_into_nuc,   'dNuc'], df.loc[sig_into_nuc,   'dChrom'],
               s=22, color=C_NUC,   alpha=0.80, linewidths=0, rasterized=True, zorder=4,
               label=f'INTO Nuclear  (n = {n_nuc})')
    ax.scatter(df.loc[sig_into_chrom, 'dNuc'], df.loc[sig_into_chrom, 'dChrom'],
               s=22, color=C_CHROM, alpha=0.80, linewidths=0, rasterized=True, zorder=5,
               label=f'INTO Chromatin  (n = {n_chrom})')

    # Diagonal region labels — placed well inside axes bounds
    ax.text(0.08, 0.78, 'INTO\nChromatin', transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=C_CHROM,
            ha='center', va='center', alpha=0.85)
    ax.text(0.88, 0.20, 'INTO\nNuclear', transform=ax.transAxes,
            fontsize=11, fontweight='bold', color=C_NUC,
            ha='center', va='center', alpha=0.85)

    # Correlation stats box — bottom-left, away from dense dots
    stat_txt = (f'Pearson r = {r_p:.3f}  (p = {p_p:.1e})\n'
                f'Spearman ρ = {r_s:.3f}  (p = {p_s:.1e})')
    ax.text(0.03, 0.03, stat_txt, transform=ax.transAxes,
            fontsize=9.5, va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                      edgecolor='#BBBBBB', alpha=0.95, zorder=10),
            zorder=10)

    # Axes
    ax.set_xlabel('ΔSoluble Nuclear  (AW − Naïve)',  fontsize=12, labelpad=8)
    ax.set_ylabel('ΔChromatin  (AW − Naïve)',         fontsize=12, labelpad=8)
    ax.set_title('Chromatin vs Soluble Nuclear proteome changes\nduring Acute Withdrawal',
                 fontsize=13, fontweight='bold', pad=12)
    ax.tick_params(labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=10, framealpha=0.95, loc='upper right',
              edgecolor='#CCCCCC', handletextpad=0.5, borderpad=0.7,
              bbox_to_anchor=(0.99, 0.99))

    plt.tight_layout()
    pdf.savefig(fig, dpi=200, bbox_inches='tight')
    plt.close(fig)

print(f'Saved: {OUT_PDF}')
