import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

FILE      = '../EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_XLSX  = 'Translocation_AW_stats.xlsx'
OUT_PDF   = 'Translocation_AW_plots.pdf'

CORR_THRESH = 3.3
FC_THRESH   = 0.5

COND_SUFFIXES = {
    'Naive':                 ['N-F-1','N-F-2','N-F-3','N-M-1','N-M-2'],
    'Intoxication':          ['I-F-1','I-F-2','I-F-3','I-M-1','I-M-2'],
    'Acute Withdrawal':      ['AW-F-1','AW-F-2','AW-M-1','AW-M-2','AW-M-3'],
    'Protracted Abstinence': ['PA-F-1','PA-F-2','PA-F-3','PA-M-1','PA-M-2'],
}
COND_ORDER  = ['Naive', 'Intoxication', 'Acute Withdrawal', 'Protracted Abstinence']
AW_SUFFIXES = COND_SUFFIXES['Acute Withdrawal']

def calc_aw_stats(df, naive_cols, aw_cols):
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
    naive_cols = [c for c in df.columns if prefix + 'N-' in str(c)]
    aw_cols    = [c for c in df.columns if any(prefix + s in str(c) for s in AW_SUFFIXES)]
    cond_cols_all = {cname: [c for c in df.columns
                              if any(prefix + s in str(c) for s in suffixes)]
                     for cname, suffixes in COND_SUFFIXES.items()}
    return df, naive_cols, aw_cols, cond_cols_all

def get_values(lookup, gene, cols):
    if gene not in lookup.index:
        return np.full(len(cols), np.nan)
    row = lookup.loc[gene]
    return np.array([float(pd.to_numeric(row[c], errors='coerce'))
                     if c in row.index else np.nan for c in cols])

# ── Load ──
nuc_df,   nuc_naive,   nuc_aw,   nuc_cc   = load_compartment('Soluble nuclear', 'Nuc_')
chrom_df, chrom_naive, chrom_aw, chrom_cc = load_compartment('Chromatin',       'Chrom_')

nuc_fc,   nuc_corr   = calc_aw_stats(nuc_df,   nuc_naive,   nuc_aw)
chrom_fc, chrom_corr = calc_aw_stats(chrom_df, chrom_naive, chrom_aw)

# Sig masks
nuc_sig_mask   = (nuc_corr   > CORR_THRESH) & (nuc_fc.abs()   > FC_THRESH)
chrom_sig_mask = (chrom_corr > CORR_THRESH) & (chrom_fc.abs() > FC_THRESH)
nuc_sig_genes   = set(nuc_df.loc[nuc_sig_mask,    'Gene symbol'].astype(str))
chrom_sig_genes = set(chrom_df.loc[chrom_sig_mask, 'Gene symbol'].astype(str))
union_genes     = sorted(nuc_sig_genes | chrom_sig_genes)

# FC/corrected lookup per gene for both compartments
nuc_fc_lookup   = dict(zip(nuc_df['Gene symbol'].astype(str),   nuc_fc))
nuc_corr_lookup = dict(zip(nuc_df['Gene symbol'].astype(str),   nuc_corr))
chrom_fc_lookup = dict(zip(chrom_df['Gene symbol'].astype(str), chrom_fc))
chrom_corr_lookup = dict(zip(chrom_df['Gene symbol'].astype(str), chrom_corr))

nuc_lookup   = (nuc_df.drop_duplicates('Gene symbol')
                .set_index(nuc_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))
chrom_lookup = (chrom_df.drop_duplicates('Gene symbol')
                .set_index(chrom_df.drop_duplicates('Gene symbol')['Gene symbol'].astype(str)))

# ── Compute interaction scores ──
rows = []
for gene in union_genes:
    nuc_n   = get_values(nuc_lookup,   gene, nuc_cc['Naive'])
    nuc_aw_ = get_values(nuc_lookup,   gene, nuc_cc['Acute Withdrawal'])
    chr_n   = get_values(chrom_lookup, gene, chrom_cc['Naive'])
    chr_aw  = get_values(chrom_lookup, gene, chrom_cc['Acute Withdrawal'])

    delta_nuc   = np.nanmean(nuc_aw_) - np.nanmean(nuc_n)
    delta_chrom = np.nanmean(chr_aw)  - np.nanmean(chr_n)
    interaction = delta_chrom - delta_nuc

    # Per-protein p-value
    chrom_diffs = chr_aw[np.isfinite(chr_aw)]   - np.nanmean(chr_n)
    nuc_diffs   = nuc_aw_[np.isfinite(nuc_aw_)] - np.nanmean(nuc_n)
    if len(chrom_diffs) >= 2 and len(nuc_diffs) >= 2:
        _, p = stats.ttest_ind(chrom_diffs, nuc_diffs, equal_var=False)
    else:
        p = np.nan

    rows.append({
        'Gene':               gene,
        'In_Nuc_AW_sig':      gene in nuc_sig_genes,
        'In_Chrom_AW_sig':    gene in chrom_sig_genes,
        # Chromatin AW stats
        'Chrom_AW_FC':        round(chrom_fc_lookup.get(gene, np.nan), 4),
        'Chrom_AW_corrected': round(chrom_corr_lookup.get(gene, np.nan), 4),
        # Soluble Nuclear AW stats
        'Nuc_AW_FC':          round(nuc_fc_lookup.get(gene, np.nan), 4),
        'Nuc_AW_corrected':   round(nuc_corr_lookup.get(gene, np.nan), 4),
        # Interaction
        'delta_Chrom_AW':     round(delta_chrom, 4) if np.isfinite(delta_chrom) else np.nan,
        'delta_Nuc_AW':       round(delta_nuc,   4) if np.isfinite(delta_nuc)   else np.nan,
        'Interaction_AW':     round(interaction, 4) if np.isfinite(interaction)  else np.nan,
        'pval_raw':           p,
    })

df_stats = pd.DataFrame(rows)

# BH correction
valid = df_stats['pval_raw'].notna()
pv    = df_stats.loc[valid, 'pval_raw'].values
ranks = pd.Series(pv).rank(ascending=True).values
bh    = np.minimum(1, pv * len(pv) / ranks)
for i in range(len(bh) - 2, -1, -1):
    bh[i] = min(bh[i], bh[i + 1])
df_stats.loc[valid, 'pval_BH'] = bh

df_stats = df_stats.sort_values('Interaction_AW', ascending=False).reset_index(drop=True)

# ── Sheet 1: 429 proteins shifted INTO Chromatin (for Metascape) ──
into_chrom = df_stats[(df_stats['Interaction_AW'] > 0) & (df_stats['pval_BH'] < 0.05)].copy()
into_nuc   = df_stats[(df_stats['Interaction_AW'] < 0) & (df_stats['pval_BH'] < 0.05)].copy()
print(f'Into Chromatin (BH<0.05): {len(into_chrom)}')
print(f'Into Nuclear   (BH<0.05): {len(into_nuc)}')

# ── Cluster tests (non-circular) ──
# Cluster A: sig UP in Chromatin AW → test if their Nuclear ALSO went down
chrom_up_genes = set(chrom_df.loc[
    (chrom_corr > CORR_THRESH) & (chrom_fc > FC_THRESH), 'Gene symbol'].astype(str))
clust_A_nuc_deltas = []
for gene in chrom_up_genes:
    nuc_n   = get_values(nuc_lookup, gene, nuc_cc['Naive'])
    nuc_aw_ = get_values(nuc_lookup, gene, nuc_cc['Acute Withdrawal'])
    d = np.nanmean(nuc_aw_) - np.nanmean(nuc_n)
    if np.isfinite(d):
        clust_A_nuc_deltas.append(d)
clust_A_nuc_deltas = np.array(clust_A_nuc_deltas)
t_A, p_A = stats.ttest_1samp(clust_A_nuc_deltas, popmean=0)
print(f'\nCluster A (Chrom UP, n={len(chrom_up_genes)}):')
print(f'  Mean ΔNuclear = {np.mean(clust_A_nuc_deltas):.4f}  t={t_A:.3f}  p={p_A:.4e}')
print(f'  → {"Nuclear DECREASED" if t_A < 0 else "Nuclear INCREASED"} in these proteins (p={p_A:.4e})')

# Cluster B: sig DOWN in Nuclear AW → test if their Chromatin went up
nuc_down_genes = set(nuc_df.loc[
    (nuc_corr > CORR_THRESH) & (nuc_fc < -FC_THRESH), 'Gene symbol'].astype(str))
clust_B_chrom_deltas = []
for gene in nuc_down_genes:
    chr_n  = get_values(chrom_lookup, gene, chrom_cc['Naive'])
    chr_aw = get_values(chrom_lookup, gene, chrom_cc['Acute Withdrawal'])
    d = np.nanmean(chr_aw) - np.nanmean(chr_n)
    if np.isfinite(d):
        clust_B_chrom_deltas.append(d)
clust_B_chrom_deltas = np.array(clust_B_chrom_deltas)
t_B, p_B = stats.ttest_1samp(clust_B_chrom_deltas, popmean=0)
print(f'\nCluster B (Nuclear DOWN, n={len(nuc_down_genes)}):')
print(f'  Mean ΔChromatin = {np.mean(clust_B_chrom_deltas):.4f}  t={t_B:.3f}  p={p_B:.4e}')
print(f'  → {"Chromatin INCREASED" if t_B > 0 else "Chromatin DECREASED"} in these proteins (p={p_B:.4e})')

# Cluster C: both — UP in Chrom AND DOWN in Nuc (strict bidirectional)
strict_transloc = chrom_up_genes & set(nuc_df.loc[
    (nuc_corr > CORR_THRESH) & (nuc_fc < -FC_THRESH), 'Gene symbol'].astype(str))
print(f'\nStrict bidirectional (UP in Chrom AND DOWN in Nuc): n={len(strict_transloc)}')

# ── Save Excel ──
cluster_summary = pd.DataFrame([
    {'Test': 'Cluster A: Chrom-UP genes — is mean ΔNuclear < 0?',
     'N': len(chrom_up_genes), 'Mean_delta': round(np.mean(clust_A_nuc_deltas), 4),
     't_stat': round(t_A, 4), 'p_value': p_A,
     'Interpretation': 'Nuclear decreased' if t_A < 0 else 'Nuclear increased'},
    {'Test': 'Cluster B: Nuc-DOWN genes — is mean ΔChromatin > 0?',
     'N': len(nuc_down_genes), 'Mean_delta': round(np.mean(clust_B_chrom_deltas), 4),
     't_stat': round(t_B, 4), 'p_value': p_B,
     'Interpretation': 'Chromatin increased' if t_B > 0 else 'Chromatin decreased'},
])

with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as writer:
    into_chrom.to_excel(writer, sheet_name='Into_Chromatin_AW', index=False)
    into_nuc.to_excel(writer,   sheet_name='Into_Nuclear_AW',   index=False)
    df_stats.to_excel(writer,   sheet_name='All_proteins',       index=False)
    cluster_summary.to_excel(writer, sheet_name='Cluster_tests', index=False)
    pd.DataFrame(sorted(strict_transloc), columns=['Gene']).to_excel(
        writer, sheet_name='Strict_bidirectional', index=False)
print(f'Saved: {OUT_XLSX}')

# ── Plots ──
with PdfPages(OUT_PDF) as pdf:

    # 1. Scatter: delta_Chrom vs delta_Nuc, colored by interaction
    fig, ax = plt.subplots(figsize=(8, 7))
    fig.patch.set_facecolor('white')
    ia = df_stats['Interaction_AW'].values
    dc = df_stats['delta_Chrom_AW'].values
    dn = df_stats['delta_Nuc_AW'].values
    fin = np.isfinite(ia) & np.isfinite(dc) & np.isfinite(dn)
    sc = ax.scatter(dn[fin], dc[fin], c=ia[fin], cmap='RdBu_r', vmin=-2, vmax=2,
                    s=14, alpha=0.6, linewidths=0, rasterized=True)
    ax.axline((0,0), slope=1, color='#888888', linestyle='--', linewidth=1.2, label='Chrom = Nuc')
    ax.axhline(0, color='#cccccc', linewidth=0.8)
    ax.axvline(0, color='#cccccc', linewidth=0.8)
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label('Interaction score\n(ΔChrom − ΔNuc)', fontsize=10)
    ax.set_xlabel('ΔNuclear  (AW − Naïve)', fontsize=12)
    ax.set_ylabel('ΔChromatin  (AW − Naïve)', fontsize=12)
    ax.set_title('Protein translocation: Chromatin vs Nuclear during AW\n'
                 'Points above diagonal → shifted INTO Chromatin', fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(fontsize=10)
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight'); plt.close(fig)

    # 2. Volcano: interaction score vs -log10(p_raw)
    fig, ax = plt.subplots(figsize=(8, 7))
    fig.patch.set_facecolor('white')
    ia2  = df_stats['Interaction_AW'].values
    praw = df_stats['pval_raw'].values
    pbh  = df_stats['pval_BH'].values
    logp = -np.log10(np.where(praw > 0, praw, np.nan))
    fin2 = np.isfinite(ia2) & np.isfinite(logp)
    su = fin2 & (ia2 > 0) & (pbh < 0.05)
    sd = fin2 & (ia2 < 0) & (pbh < 0.05)
    ns = fin2 & ~su & ~sd
    ax.scatter(ia2[ns],  logp[ns],  s=12, color='#CCCCCC', alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(ia2[sd],  logp[sd],  s=25, color='#2B7FD4', alpha=0.85, linewidths=0, label=f'Into Nuclear (n={sd.sum()})')
    ax.scatter(ia2[su],  logp[su],  s=25, color='#E8305A', alpha=0.85, linewidths=0, label=f'Into Chromatin (n={su.sum()})')
    ax.axhline(-np.log10(0.05), color='#555555', linestyle='--', linewidth=1.2, alpha=0.7)
    ax.axvline(0, color='#999999', linewidth=0.8)
    ax.set_xlabel('Interaction score  (ΔChromatin − ΔNuclear during AW)', fontsize=12)
    ax.set_ylabel('−log₁₀(p-value)', fontsize=12)
    ax.set_title('Translocation volcano — Acute Withdrawal', fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.legend(fontsize=10, framealpha=0.85)
    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight'); plt.close(fig)

    # 3. Cluster tests bar chart
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.patch.set_facecolor('white')

    # Cluster A: distribution of ΔNuclear for Chrom-UP proteins
    ax = axes[0]
    ax.hist(clust_A_nuc_deltas, bins=40, color='#E8305A', alpha=0.7, edgecolor='white')
    ax.axvline(0, color='#333333', linewidth=1.5, linestyle='--')
    ax.axvline(np.mean(clust_A_nuc_deltas), color='#B2182B', linewidth=2.5, label=f'Mean={np.mean(clust_A_nuc_deltas):.3f}')
    p_str = f'p = {p_A:.3e}' if p_A >= 1e-4 else f'p = {p_A:.2e}'
    ax.text(0.97, 0.97, p_str, transform=ax.transAxes, ha='right', va='top',
            fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#CCCCCC'))
    ax.set_xlabel('ΔNuclear  (AW − Naïve)', fontsize=11)
    ax.set_ylabel('Number of proteins', fontsize=11)
    ax.set_title(f'Cluster A: proteins UP in Chromatin (n={len(chrom_up_genes)})\n'
                 f'Did their Nuclear levels decrease during AW?', fontsize=10, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    # Cluster B: distribution of ΔChromatin for Nuc-DOWN proteins
    ax = axes[1]
    ax.hist(clust_B_chrom_deltas, bins=40, color='#2B7FD4', alpha=0.7, edgecolor='white')
    ax.axvline(0, color='#333333', linewidth=1.5, linestyle='--')
    ax.axvline(np.mean(clust_B_chrom_deltas), color='#2166AC', linewidth=2.5, label=f'Mean={np.mean(clust_B_chrom_deltas):.3f}')
    p_str = f'p = {p_B:.3e}' if p_B >= 1e-4 else f'p = {p_B:.2e}'
    ax.text(0.97, 0.97, p_str, transform=ax.transAxes, ha='right', va='top',
            fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#CCCCCC'))
    ax.set_xlabel('ΔChromatin  (AW − Naïve)', fontsize=11)
    ax.set_ylabel('Number of proteins', fontsize=11)
    ax.set_title(f'Cluster B: proteins DOWN in Nuclear (n={len(nuc_down_genes)})\n'
                 f'Did their Chromatin levels increase during AW?', fontsize=10, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    plt.tight_layout()
    pdf.savefig(fig, dpi=150, bbox_inches='tight'); plt.close(fig)

print(f'Saved: {OUT_PDF}')
