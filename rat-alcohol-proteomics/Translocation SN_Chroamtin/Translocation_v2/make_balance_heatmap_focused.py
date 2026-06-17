import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

FILE = 'Chrom_vs_SN_Balance.xlsx'
OUT  = 'Balance_Heatmap_Focused.pdf'

df = pd.read_excel(FILE, sheet_name='Chrom_vs_SN_Balance')
both_sig = df[(df['Sig_Naive'] != 'NS') & (df['Sig_AW'] != 'NS')].copy()

CMAP = 'RdBu_r'

def draw(ax, mat, genes, title, norm, fontsize=7):
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Naive\n(Chrom−SN)', 'AW\n(Chrom−SN)'], fontsize=9, fontweight='bold')
    ax.set_title(title, fontsize=9, fontweight='bold', pad=6)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=fontsize)
    return im

with PdfPages(OUT) as pdf:

    # ── VERSION 1: Top 50 by |FC_AW|, sorted descending ─────────────────────
    top50 = both_sig.copy()
    top50['abs_AW'] = top50['FC_AW'].abs()
    top50 = top50.nlargest(50, 'abs_AW').sort_values('FC_AW', ascending=False)
    mat50 = top50[['FC_Naive', 'FC_AW']].values.astype(float)
    genes50 = top50['Gene symbol'].fillna(top50['Accession']).tolist()
    vmax50 = min(np.nanpercentile(np.abs(mat50), 99), 10.0)
    norm50 = TwoSlopeNorm(vmin=-vmax50, vcenter=0, vmax=vmax50)

    fig, axes = plt.subplots(1, 2, figsize=(7, 10),
                             gridspec_kw={'width_ratios': [4, 0.2]})
    im = draw(axes[0], mat50, genes50,
              f'Top 50 proteins by |AW FC|\n(significant in both Naive & AW,  AW-M-3 excluded)',
              norm50, fontsize=7)
    plt.colorbar(im, cax=axes[1], label='FC (Chrom − SN, log2)')
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    print(f'Version 1 (top 50): done')

    # ── VERSION 2: Switchers — opposite direction Naive vs AW ────────────────
    switchers = both_sig[
        ((both_sig['Sig_Naive'] == 'Up')   & (both_sig['Sig_AW'] == 'Down')) |
        ((both_sig['Sig_Naive'] == 'Down') & (both_sig['Sig_AW'] == 'Up'))
    ].copy()
    print(f'Switchers (direction flip Naive→AW): {len(switchers)}')
    print(f'  Naive Up → AW Down (SN→Chrom shift): {((switchers["Sig_Naive"]=="Up") & (switchers["Sig_AW"]=="Down")).sum()}')
    print(f'  Naive Down → AW Up (Chrom→SN shift): {((switchers["Sig_Naive"]=="Down") & (switchers["Sig_AW"]=="Up")).sum()}')

    # Sort: Naive-Up→AW-Down first (SN-to-Chrom), then Naive-Down→AW-Up
    sw_sn_to_chrom = switchers[(switchers['Sig_Naive']=='Up') & (switchers['Sig_AW']=='Down')].sort_values('FC_AW', ascending=False)
    sw_chrom_to_sn = switchers[(switchers['Sig_Naive']=='Down') & (switchers['Sig_AW']=='Up')].sort_values('FC_AW', ascending=False)
    switchers_sorted = pd.concat([sw_sn_to_chrom, sw_chrom_to_sn])

    mat_sw = switchers_sorted[['FC_Naive', 'FC_AW']].values.astype(float)
    genes_sw = switchers_sorted['Gene symbol'].fillna(switchers_sorted['Accession']).tolist()
    vmax_sw = min(np.nanpercentile(np.abs(mat_sw), 99), 8.0)
    norm_sw = TwoSlopeNorm(vmin=-vmax_sw, vcenter=0, vmax=vmax_sw)

    fontsize_sw = max(4, min(7, int(200 / max(len(genes_sw), 1))))
    fig_h = max(6, len(genes_sw) * 0.16 + 2)

    fig, axes = plt.subplots(1, 2, figsize=(7, fig_h),
                             gridspec_kw={'width_ratios': [4, 0.2]})
    im = draw(axes[0], mat_sw, genes_sw,
              f'Switchers: direction flip between Naive and AW\n'
              f'(n={len(switchers)},  AW-M-3 excluded)\n'
              f'Top = Naive Chrom-enriched → AW SN-enriched  |  Bottom = opposite',
              norm_sw, fontsize=fontsize_sw)
    # Divider line between the two groups
    n_first = len(sw_sn_to_chrom)
    axes[0].axhline(y=n_first - 0.5, color='white', linewidth=2, linestyle='--')
    # Group labels
    axes[0].text(1.15, (n_first/2) / len(genes_sw), f'Naive Chrom>\nAW SN>\n(n={n_first})',
                 transform=axes[0].transAxes, fontsize=6, va='center', color='steelblue')
    axes[0].text(1.15, (n_first + len(sw_chrom_to_sn)/2) / len(genes_sw),
                 f'Naive SN>\nAW Chrom>\n(n={len(sw_chrom_to_sn)})',
                 transform=axes[0].transAxes, fontsize=6, va='center', color='firebrick')

    plt.colorbar(im, cax=axes[1], label='FC (Chrom − SN, log2)')
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    print(f'Version 2 (switchers): done')

print(f'\nSaved: {OUT}')
