"""
IPA Upstream Regulator bar plot — Acute Withdrawal
===================================================
Top 10 IPA-predicted upstream regulators, ranked by -log10(p-value).
Colors: pastel red = activated (z>0), pastel blue = inhibited (z<0), grey = z=NA.

Output: IPA_AW_Upstream_Regulators.pdf / .png
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PDF = os.path.join(SCRIPT_DIR, 'IPA_AW_Upstream_Regulators.pdf')
OUT_PNG = os.path.join(SCRIPT_DIR, 'IPA_AW_Upstream_Regulators.png')

# Values read from figure — (gene, -log10 p-value, z-score or None, concordant in AW)
# Concordant = IPA direction matches proteomics direction in AW
DATA = [
    ('KRAS',    8.52, -2.71, True),
    ('MAPT',    6.20, -1.41, True),
    ('RRAS2',   5.88,  2.67, False),
    ('GLB1',    4.68, -2.24, True),
    ('SOD1',    4.22,  None, False),
    ('MRTFB',   3.52,  0.44, False),
    ('ELAVL1',  3.38, -1.61, False),
    ('CREB1',   3.05, -1.36, False),
    ('MAP2K1',  2.82,  0.71, False),
    ('HSP90B1', 2.58,  0.28, False),
]

# Pastel colors
C_ACT  = '#E8AAAA'   # pastel red  — activated (z > 0)
C_INH  = '#A8BED8'   # pastel blue — inhibited (z < 0)
C_NA   = '#C8C8C8'   # grey        — z = NA
C_EDGE = '#555555'

def bar_color(z):
    if z is None:  return C_NA
    if z > 0:      return C_ACT
    return C_INH

genes       = [d[0] for d in DATA]
pvals       = [d[1] for d in DATA]
zscores     = [d[2] for d in DATA]
concordant  = [d[3] for d in DATA]
colors      = [bar_color(z) for z in zscores]
z_labels    = [f'z={z:.2f}' if z is not None else 'z=NA' for z in zscores]
# Add asterisk to gene label if concordant
gene_labels = [f'* {g}' if c else g for g, c in zip(genes, concordant)]

fig, ax = plt.subplots(figsize=(26, 14))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

y_pos = np.arange(len(genes))

bars = ax.barh(y_pos, pvals, color=colors, edgecolor=C_EDGE,
               linewidth=1.2, height=0.65)

# Z-score labels inside bars
for i, (bar, zlbl) in enumerate(zip(bars, z_labels)):
    x_end = bar.get_width()
    ax.text(x_end - 0.20, i, zlbl,
            ha='right', va='center',
            fontsize=40, fontfamily='Arial',
            fontweight='bold', color='#333333')

ax.set_yticks(y_pos)
ax.set_yticklabels(gene_labels, fontsize=44, fontfamily='Arial', fontstyle='italic')
ax.set_xlabel('−log$_{10}$ (IPA p-value)', fontsize=44, fontfamily='Arial', labelpad=10)
ax.set_title('Upstream Regulators (RNA-seq × Proteomics)', fontsize=48, fontweight='bold',
             fontfamily='Arial', pad=20)

ax.set_xlim(0, max(pvals) * 1.08)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2.0)
ax.spines['bottom'].set_linewidth(2.0)
ax.tick_params(axis='x', labelsize=40)
ax.tick_params(axis='y', length=0, pad=8)

# Legend — inside plot, bottom right, includes asterisk explanation
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor=C_ACT, edgecolor=C_EDGE, linewidth=1.2, label='Activated (z > 0)'),
    Patch(facecolor=C_INH, edgecolor=C_EDGE, linewidth=1.2, label='Inhibited (z < 0)'),
    Patch(facecolor=C_NA,  edgecolor=C_EDGE, linewidth=1.2, label='z = NA'),
    Line2D([0], [0], marker='', color='none',
           label='* Concordant with AW proteomics'),
]
leg = ax.legend(handles=legend_elements, fontsize=36, framealpha=0.95,
                edgecolor='#CCCCCC', loc='lower right', handlelength=1.2,
                borderpad=0.8, labelspacing=0.5)

plt.subplots_adjust(left=0.18, right=0.97, top=0.92, bottom=0.12)

with PdfPages(OUT_PDF) as pdf:
    pdf.savefig(fig, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {OUT_PDF}')
print(f'Saved: {OUT_PNG}')
