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
OUT  = 'Balance_Heatmap_Naive_vs_AW.pdf'

df = pd.read_excel(FILE, sheet_name='Chrom_vs_SN_Balance')

# Keep only proteins significant in BOTH Naive and AW
both_sig = df[(df['Sig_Naive'] != 'NS') & (df['Sig_AW'] != 'NS')].copy()
print(f'Significant in both Naive and AW: {len(both_sig)}')
print(f'  Naive Up / Down: {(both_sig["Sig_Naive"]=="Up").sum()} / {(both_sig["Sig_Naive"]=="Down").sum()}')
print(f'  AW   Up / Down: {(both_sig["Sig_AW"]=="Up").sum()} / {(both_sig["Sig_AW"]=="Down").sum()}')

mat = both_sig[['FC_Naive', 'FC_AW']].values.astype(float)
genes = both_sig['Gene symbol'].fillna(both_sig['Accession']).tolist()

# Cluster rows
if len(both_sig) >= 3:
    dist = pdist(mat, metric='euclidean')
    link = linkage(dist, method='ward')
    order = leaves_list(link)
else:
    order = np.arange(len(both_sig))

mat   = mat[order]
genes = [genes[i] for i in order]

vmax = min(np.nanpercentile(np.abs(mat), 98), 5.0)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

show_labels = len(genes) <= 150
fig_h = max(6, len(genes) * 0.13 + 2) if show_labels else 10
fig_w = 6

fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h),
                          gridspec_kw={'width_ratios': [4, 0.25]})
ax, cax = axes

im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', norm=norm, interpolation='nearest')

ax.set_xticks([0, 1])
ax.set_xticklabels(['Naive\n(Chrom − SN)', 'AW\n(Chrom − SN)'], fontsize=9, fontweight='bold')
ax.set_title(
    f'Chrom vs SN balance: Naive vs Acute Withdrawal\n'
    f'Significant in both conditions  (n={len(genes)})\n'
    f'AW-M-3 excluded',
    fontsize=9, fontweight='bold', pad=8)

if show_labels:
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=5)
else:
    ax.set_yticks([])
    ax.set_ylabel(f'{len(genes)} proteins (labels hidden — too many to display)', fontsize=7)

plt.colorbar(im, cax=cax, label='FC (Chrom − SN, log2)')
cax.tick_params(labelsize=7)

plt.tight_layout()
with PdfPages(OUT) as pdf:
    pdf.savefig(fig, bbox_inches='tight')
plt.close()

print(f'Saved: {OUT}')
