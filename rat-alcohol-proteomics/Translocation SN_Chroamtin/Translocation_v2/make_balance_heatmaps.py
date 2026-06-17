import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
from openpyxl import load_workbook

FILE = 'Chrom_vs_SN_Balance.xlsx'
OUT  = 'Balance_Heatmaps.pdf'

CONDS = ['Naive', 'Intox', 'AW', 'PA']
FC_COLS   = [f'FC_{c}'        for c in CONDS]
SIG_COLS  = [f'Sig_{c}'       for c in CONDS]
CORR_COLS = [f'Corrected_{c}' for c in CONDS]

df = pd.read_excel(FILE, sheet_name='Chrom_vs_SN_Balance')

# Sig helpers
def sig_in(row, conds):
    return any(row[f'Sig_{c}'] != 'NS' for c in conds)

# Subsets
aw_sig   = df[df['Sig_AW'] != 'NS'].copy()
any_sig  = df[df[SIG_COLS].apply(lambda r: any(r != 'NS'), axis=1)].copy()
aw_up    = df[df['Sig_AW'] == 'Up'].copy()    # more in Chrom during AW
aw_down  = df[df['Sig_AW'] == 'Down'].copy()  # more in SN during AW

CMAP = 'RdBu_r'   # red = Chrom-enriched, blue = SN-enriched

def fc_matrix(subset):
    mat = subset[FC_COLS].values.astype(float)
    mat = np.nan_to_num(mat, nan=0.0)
    return mat

def make_norm(mat, vcap=3.0):
    vmax = min(np.nanpercentile(np.abs(mat), 98), vcap)
    return TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

def draw_heatmap(ax, mat, genes, title, norm, show_genes=True, fontsize=5):
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax.set_xticks(range(len(CONDS)))
    ax.set_xticklabels(CONDS, fontsize=8, fontweight='bold')
    ax.set_title(title, fontsize=9, fontweight='bold', pad=6)
    if show_genes and len(genes) <= 120:
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=fontsize)
    else:
        ax.set_yticks([])
    ax.set_ylabel(f'n = {len(genes)} proteins', fontsize=7)
    return im

def cluster_rows(mat):
    if len(mat) < 3:
        return np.arange(len(mat))
    dist = pdist(mat, metric='euclidean')
    link = linkage(dist, method='ward')
    return leaves_list(link)

from matplotlib.backends.backend_pdf import PdfPages

with PdfPages(OUT) as pdf:

    # ── PAGE 1: AW-significant, sorted by AW FC ──────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, max(6, len(aw_sig)*0.13 + 2)),
                             gridspec_kw={'width_ratios': [4, 0.15]})
    ax, cax = axes
    sub = aw_sig.sort_values('FC_AW', ascending=False)
    mat = fc_matrix(sub)
    norm = make_norm(mat)
    im = draw_heatmap(ax, mat, sub['Gene symbol'].fillna(sub['Accession']).tolist(),
                      f'AW-significant proteins sorted by AW FC\n'
                      f'(Chrom − SN, positive = Chrom-enriched)  n={len(sub)}',
                      norm, fontsize=5)
    plt.colorbar(im, cax=cax, label='FC (Chrom − SN)')
    fig.suptitle('Option 1: AW-significant, sorted by AW FC', fontsize=11, fontweight='bold', y=1.01)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # ── PAGE 2: AW-significant, clustered ────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, max(6, len(aw_sig)*0.13 + 2)),
                             gridspec_kw={'width_ratios': [4, 0.15]})
    ax, cax = axes
    mat = fc_matrix(aw_sig)
    order = cluster_rows(mat)
    mat_c = mat[order]
    genes_c = aw_sig['Gene symbol'].fillna(aw_sig['Accession']).iloc[order].tolist()
    norm = make_norm(mat_c)
    im = draw_heatmap(ax, mat_c, genes_c,
                      f'AW-significant proteins — hierarchical clustering\n'
                      f'(Chrom − SN)  n={len(aw_sig)}',
                      norm, fontsize=5)
    plt.colorbar(im, cax=cax, label='FC (Chrom − SN)')
    fig.suptitle('Option 2: AW-significant, clustered', fontsize=11, fontweight='bold', y=1.01)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # ── PAGE 3: AW-sig split Up / Down side by side ──────────────────────────
    n_up, n_dn = len(aw_up), len(aw_down)
    fig_h = max(8, max(n_up, n_dn)*0.13 + 3)
    fig, axes = plt.subplots(1, 4, figsize=(14, fig_h),
                             gridspec_kw={'width_ratios': [4, 0.15, 4, 0.15]})
    ax1, cax1, ax2, cax2 = axes

    sub_up = aw_up.sort_values('FC_AW', ascending=False)
    mat_up = fc_matrix(sub_up)
    norm_up = make_norm(mat_up)
    im1 = draw_heatmap(ax1, mat_up,
                       sub_up['Gene symbol'].fillna(sub_up['Accession']).tolist(),
                       f'AW Up (Chrom > SN)  n={n_up}', norm_up, fontsize=5)
    plt.colorbar(im1, cax=cax1, label='FC')

    sub_dn = aw_down.sort_values('FC_AW', ascending=True)
    mat_dn = fc_matrix(sub_dn)
    norm_dn = make_norm(mat_dn)
    im2 = draw_heatmap(ax2, mat_dn,
                       sub_dn['Gene symbol'].fillna(sub_dn['Accession']).tolist(),
                       f'AW Down (SN > Chrom)  n={n_dn}', norm_dn, fontsize=5)
    plt.colorbar(im2, cax=cax2, label='FC')

    fig.suptitle('Option 3: AW-significant split — Up vs Down', fontsize=11, fontweight='bold')
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # ── PAGE 4: Any-condition significant, clustered ──────────────────────────
    sub_any = any_sig.copy()
    mat_any = fc_matrix(sub_any)
    order_any = cluster_rows(mat_any)
    mat_any_c = mat_any[order_any]
    show_genes = len(sub_any) <= 150
    fig_h = max(8, len(sub_any)*0.10 + 2) if show_genes else 10
    fig, axes = plt.subplots(1, 2, figsize=(12, fig_h),
                             gridspec_kw={'width_ratios': [4, 0.15]})
    ax, cax = axes
    genes_any = sub_any['Gene symbol'].fillna(sub_any['Accession']).iloc[order_any].tolist()
    norm_any = make_norm(mat_any_c)
    im = draw_heatmap(ax, mat_any_c, genes_any,
                      f'All significant in any condition — clustered\n'
                      f'(Chrom − SN)  n={len(sub_any)}',
                      norm_any, show_genes=show_genes, fontsize=4)
    plt.colorbar(im, cax=cax, label='FC (Chrom − SN)')
    fig.suptitle('Option 4: Any-condition significant, clustered', fontsize=11, fontweight='bold', y=1.01)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # ── PAGE 5: AW-significant only, show corrected-p as separate panel ──────
    sub = aw_sig.sort_values('FC_AW', ascending=False)
    mat_fc   = fc_matrix(sub)
    mat_corr = sub[CORR_COLS].values.astype(float)
    mat_corr = np.nan_to_num(mat_corr, nan=0.0)

    fig_h = max(6, len(sub)*0.13 + 2)
    fig, axes = plt.subplots(1, 4, figsize=(16, fig_h),
                             gridspec_kw={'width_ratios': [4, 0.15, 4, 0.15]})
    ax1, cax1, ax2, cax2 = axes

    norm_fc = make_norm(mat_fc)
    im1 = draw_heatmap(ax1, mat_fc,
                       sub['Gene symbol'].fillna(sub['Accession']).tolist(),
                       f'FC (Chrom − SN)  n={len(sub)}', norm_fc, fontsize=5)
    plt.colorbar(im1, cax=cax1, label='FC')

    from matplotlib.colors import Normalize
    norm_p = Normalize(vmin=0, vmax=np.nanpercentile(mat_corr, 99))
    im2 = ax2.imshow(mat_corr, aspect='auto', cmap='YlOrRd', norm=norm_p, interpolation='nearest')
    ax2.set_xticks(range(len(CONDS)))
    ax2.set_xticklabels(CONDS, fontsize=8, fontweight='bold')
    ax2.set_title(f'Corrected −log2(p)  n={len(sub)}', fontsize=9, fontweight='bold', pad=6)
    ax2.set_yticks([])
    plt.colorbar(im2, cax=cax2, label='−log2(adj p)')
    # add threshold line indicator
    ax2.axvline(x=-0.5, color='white', linewidth=0.5)

    fig.suptitle('Option 5: AW-significant — FC and corrected-p side by side', fontsize=11, fontweight='bold', y=1.01)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

print(f'Saved: {OUT}')
print(f'\nSummary:')
print(f'  AW-significant total : {len(aw_sig)}  (Up={len(aw_up)}, Down={len(aw_down)})')
print(f'  Any-condition sig    : {len(any_sig)}')
