import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import hypergeom, pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

FILE    = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_PDF = 'RRHO2_Nuc_vs_Chrom_AW.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

COND_SUFFIXES = {
    'Naive':            ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Acute Withdrawal': ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
}

# ── Colormap: white → yellow → orange → red (classic RRHO2 style) ──
RRHO_CMAP = LinearSegmentedColormap.from_list(
    'rrho2',
    ['#FFFFFF', '#FFFACD', '#FFEDA0', '#FED976', '#FEB24C',
     '#FD8D3C', '#FC4E2A', '#E31A1C', '#BD0026', '#800026'],
    N=256
)

def calc_delta(df, naive_cols, cond_cols):
    """Mean(condition) - Mean(naive) per row."""
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    return cond.mean(axis=1) - naive.mean(axis=1)

def load_compartment(excel_name, prefix):
    raw = pd.read_excel(FILE, sheet_name=excel_name, header=0)
    df  = raw[raw['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)
    naive_cols = [c for c in df.columns
                  if any(prefix + s in str(c) for s in COND_SUFFIXES['Naive'])]
    aw_cols    = [c for c in df.columns
                  if any(prefix + s in str(c) for s in COND_SUFFIXES['Acute Withdrawal'])]
    return df, naive_cols, aw_cols

def bh_correct_matrix(pmat):
    """Benjamini-Hochberg correction across all cells of a matrix."""
    flat    = pmat.ravel()
    n       = len(flat)
    order   = np.argsort(flat)
    ranks   = np.empty_like(order); ranks[order] = np.arange(1, n+1)
    bh_flat = np.minimum(1.0, flat * n / ranks)
    # Enforce monotonicity
    sorted_bh = bh_flat[order]
    for i in range(len(sorted_bh) - 2, -1, -1):
        sorted_bh[i] = min(sorted_bh[i], sorted_bh[i+1])
    bh_flat[order] = sorted_bh
    return bh_flat.reshape(pmat.shape)

def rrho2_quadrant(genes_by_1, genes_by_2, dir1, dir2, N, steps):
    """
    Compute -log10(BH p) matrix for one RRHO2 quadrant.
    dir = 'up'   → take from front of sorted list (highest delta)
    dir = 'down' → take from back  of sorted list (lowest delta)
    """
    n   = len(steps)
    mat = np.ones((n, n))   # default p=1
    for i, s1 in enumerate(steps):
        set1 = set(genes_by_1[:s1]) if dir1 == 'up' else set(genes_by_1[-s1:])
        for j, s2 in enumerate(steps):
            set2 = set(genes_by_2[:s2]) if dir2 == 'up' else set(genes_by_2[-s2:])
            ov   = len(set1 & set2)
            p    = float(hypergeom.sf(max(0, ov - 1), N, s1, s2))
            mat[i, j] = max(p, 1e-300)
    bh = bh_correct_matrix(mat)
    return -np.log10(np.maximum(bh, 1e-300))

# ── Load data ──
nuc_df,   nuc_naive,   nuc_aw   = load_compartment('Soluble nuclear', 'Nuc_')
chrom_df, chrom_naive, chrom_aw = load_compartment('Chromatin',       'Chrom_')

delta_nuc_s   = calc_delta(nuc_df,   nuc_naive,   nuc_aw)
delta_chrom_s = calc_delta(chrom_df, chrom_naive, chrom_aw)

# Gene-level lookup
nuc_delta   = dict(zip(nuc_df['Gene symbol'].astype(str),   delta_nuc_s))
chrom_delta = dict(zip(chrom_df['Gene symbol'].astype(str), delta_chrom_s))

# Union of AW-sig proteins (reuse threshold from earlier analysis)
from scipy.stats import ttest_ind

def calc_corr(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
    fc    = cond.mean(axis=1) - naive.mean(axis=1)
    p_vals = []
    for i in range(len(df)):
        n = naive.iloc[i].dropna().values
        c = cond.iloc[i].dropna().values
        if len(n) >= 2 and len(c) >= 2:
            _, p = ttest_ind(n, c, equal_var=False)
            p_vals.append(p)
        else:
            p_vals.append(np.nan)
    p_vals  = pd.Series(p_vals, index=df.index)
    valid   = p_vals.notna()
    ranks   = p_vals[valid].rank(ascending=False)
    corr    = pd.Series(np.nan, index=df.index)
    corr[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corr

nuc_fc,   nuc_corr   = calc_corr(nuc_df,   nuc_naive,   nuc_aw)
chrom_fc, chrom_corr = calc_corr(chrom_df, chrom_naive, chrom_aw)

nuc_sig   = set(nuc_df.loc[(nuc_corr   > CORR_THRESH) & (nuc_fc.abs()   > FC_THRESH), 'Gene symbol'].astype(str))
chrom_sig = set(chrom_df.loc[(chrom_corr > CORR_THRESH) & (chrom_fc.abs() > FC_THRESH), 'Gene symbol'].astype(str))
union_genes = sorted(nuc_sig | chrom_sig)
N = len(union_genes)
print(f'Total proteins for RRHO2: {N}')

# Build delta arrays for valid genes only
valid_genes = [g for g in union_genes
               if np.isfinite(nuc_delta.get(g, np.nan))
               and np.isfinite(chrom_delta.get(g, np.nan))]
N_valid = len(valid_genes)
print(f'Proteins with both delta values: {N_valid}')

# ── Pearson/Spearman correlation (p-value for Figure 2) ──
dc = np.array([chrom_delta[g] for g in valid_genes])
dn = np.array([nuc_delta[g]   for g in valid_genes])
r_pearson,  p_pearson  = pearsonr(dc,  dn)
r_spearman, p_spearman = spearmanr(dc, dn)
print(f'\nFigure 2 p-values:')
print(f'  Pearson  r = {r_pearson:.4f},  p = {p_pearson:.4e}')
print(f'  Spearman r = {r_spearman:.4f},  p = {p_spearman:.4e}')

# ── Rank genes ──
by_chrom = sorted(valid_genes, key=lambda g: chrom_delta[g], reverse=True)
by_nuc   = sorted(valid_genes, key=lambda g: nuc_delta[g],   reverse=True)

# Steps: ~1% of N per step, targeting ~100 steps
step_size = max(1, N_valid // 100)
steps     = list(range(step_size, N_valid + 1, step_size))
n_steps   = len(steps)
print(f'RRHO2 step size: {step_size}, steps: {n_steps}')

# ── Compute 4 quadrants ──
print('Computing UU (both up)...')
log_uu = rrho2_quadrant(by_chrom, by_nuc, 'up',   'up',   N_valid, steps)
print('Computing DD (both down)...')
log_dd = rrho2_quadrant(by_chrom, by_nuc, 'down', 'down', N_valid, steps)
print('Computing UD (Chrom up, Nuc down)...')
log_ud = rrho2_quadrant(by_chrom, by_nuc, 'up',   'down', N_valid, steps)
print('Computing DU (Chrom down, Nuc up)...')
log_du = rrho2_quadrant(by_chrom, by_nuc, 'down', 'up',   N_valid, steps)

# ── Assemble full 4-quadrant matrix ──
# Layout (imshow origin='lower', y=0 at bottom):
#   Top-left    = DU (Chrom down, Nuc up)   [n:2n, 0:n]
#   Top-right   = DD (both down)             [n:2n, n:2n]
#   Bottom-left = UU (both up)               [0:n,  0:n]
#   Bottom-right= UD (Chrom up, Nuc down)    [0:n,  n:2n]
#
# Each quadrant: row/col index 0 = closest to center (smallest set)
# → need to flip so smallest sets are at center

full = np.zeros((2 * n_steps, 2 * n_steps))
full[0:n_steps,         0:n_steps]         = np.flipud(np.fliplr(log_uu))  # BL
full[0:n_steps,         n_steps:2*n_steps] = np.flipud(log_ud)             # BR
full[n_steps:2*n_steps, 0:n_steps]         = np.fliplr(log_du)             # TL
full[n_steps:2*n_steps, n_steps:2*n_steps] = log_dd                        # TR

vmax = np.nanpercentile(full[full > 0], 99)
vmax = max(vmax, 2.0)
print(f'Color scale vmax: {vmax:.2f}  (−log₁₀ BH p)')

# ── Plot ──
with PdfPages(OUT_PDF) as pdf:
    fig, ax = plt.subplots(figsize=(9, 8))
    fig.patch.set_facecolor('white')

    im = ax.imshow(full, aspect='auto', cmap=RRHO_CMAP,
                   vmin=0, vmax=vmax, origin='lower', interpolation='bilinear')

    # Quadrant dividers
    ax.axhline(n_steps - 0.5, color='#333333', linewidth=2.0)
    ax.axvline(n_steps - 0.5, color='#333333', linewidth=2.0)

    # Quadrant labels
    fs = 11
    ax.text(n_steps * 0.5,       n_steps * 0.5,       'Both UP',
            ha='center', va='center', fontsize=fs, fontweight='bold',
            color='#222222', alpha=0.45)
    ax.text(n_steps * 1.5,       n_steps * 0.5,       'Chrom UP\nNuc DOWN',
            ha='center', va='center', fontsize=fs, fontweight='bold',
            color='#222222', alpha=0.45)
    ax.text(n_steps * 0.5,       n_steps * 1.5,       'Chrom DOWN\nNuc UP',
            ha='center', va='center', fontsize=fs, fontweight='bold',
            color='#222222', alpha=0.45)
    ax.text(n_steps * 1.5,       n_steps * 1.5,       'Both DOWN',
            ha='center', va='center', fontsize=fs, fontweight='bold',
            color='#222222', alpha=0.45)

    # Axis tick labels (gene rank at edges)
    tick_pos   = [0, n_steps//2, n_steps-1, n_steps, int(n_steps*1.5), 2*n_steps-1]
    tick_labs_x = [f'{steps[-1]}\n(Nuc↑)', '', f'{steps[0]}',
                   f'{steps[0]}', '', f'{steps[-1]}\n(Nuc↓)']
    tick_labs_y = [f'{steps[-1]}\n(Chrom↑)', '', f'{steps[0]}',
                   f'{steps[0]}', '', f'{steps[-1]}\n(Chrom↓)']
    ax.set_xticks(tick_pos); ax.set_xticklabels(tick_labs_x, fontsize=8)
    ax.set_yticks(tick_pos); ax.set_yticklabels(tick_labs_y, fontsize=8)

    ax.set_xlabel('Proteins ranked by ΔSoluble Nuclear  (AW − Naïve)',
                  fontsize=12, labelpad=8)
    ax.set_ylabel('Proteins ranked by ΔChromatin  (AW − Naïve)',
                  fontsize=12, labelpad=8)
    ax.set_title(
        'RRHO2 — Chromatin vs Soluble Nuclear during Acute Withdrawal\n'
        f'Pearson r = {r_pearson:.3f}  (p = {p_pearson:.2e})   '
        f'Spearman ρ = {r_spearman:.3f}  (p = {p_spearman:.2e})',
        fontsize=12, fontweight='bold', pad=12)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label('−log₁₀(BH-corrected p)', fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close(fig)

print(f'Saved: {OUT_PDF}')
