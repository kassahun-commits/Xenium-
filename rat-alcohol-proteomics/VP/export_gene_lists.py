import pandas as pd
import numpy as np
from scipy import stats
import os

FILE = 'EDIT  Philipp Alcohol proteome Simplified - Jan 2024  copy.xlsx'
OUT_DIR = 'Metascape_lists'
os.makedirs(OUT_DIR, exist_ok=True)

CORR_THRESH = 3.3
FC_THRESH   = 0.5

def calc_stats(df, naive_cols, cond_cols):
    naive = df[naive_cols].apply(pd.to_numeric, errors='coerce')
    cond  = df[cond_cols].apply(pd.to_numeric, errors='coerce')
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
    p_vals    = pd.Series(p_vals, index=df.index)
    valid     = p_vals.notna()
    ranks     = p_vals[valid].rank(ascending=False)
    corrected = pd.Series(np.nan, index=df.index)
    corrected[valid] = -np.log2(np.minimum(1, p_vals[valid] * valid.sum() / ranks))
    return fc, corrected

CONDITIONS = [
    ('Intoxication',         'I'),
    ('Acute_Withdrawal',     'AW'),
    ('Protracted_Abstinence','PA'),
]

# (key, sheet_name_in_excel, naive_pat, cond_pat_tmpl, filter_mode)
SHEETS = [
    ('Chromatin',       'Chromatin',       'Chrom_N-', 'Chrom_{cond}-', 'keep_review'),
    ('Soluble_nuclear', 'Soluble nuclear', 'Nuc_N-',   'Nuc_{cond}-',   'keep_review'),
    ('Membrane',        'Membrane',        'Memb_N-',  'Memb_{cond}-',  'all'),
    ('Cytosol',         'Cytosol',         'Cyto_N-',  'Cyto_{cond}-',  'all'),
]

sheet_cache = {key: pd.read_excel(FILE, sheet_name=excel_name, header=0)
               for key, excel_name, _, _, _ in SHEETS}

summary = []

for sheet_key, excel_name, naive_pat, cond_pat_tmpl, fmode in SHEETS:
    for cond_label, cond_code in CONDITIONS:
        cond_pat = cond_pat_tmpl.replace('{cond}', cond_code)
        df = sheet_cache[sheet_key].copy()

        if fmode == 'keep_review':
            df = df[df['Filter'].isin(['Keep', 'Review'])].reset_index(drop=True)

        naive_cols = [c for c in df.columns if naive_pat in str(c)]
        cond_cols  = [c for c in df.columns if cond_pat  in str(c)]
        fc, corrected = calc_stats(df, naive_cols, cond_cols)

        mask = fc.notna() & corrected.notna()
        up   = mask & (corrected > CORR_THRESH) & (fc >  FC_THRESH)
        down = mask & (corrected > CORR_THRESH) & (fc < -FC_THRESH)

        genes_up   = df.loc[up,   'Gene symbol'].dropna().astype(str).tolist()
        genes_down = df.loc[down, 'Gene symbol'].dropna().astype(str).tolist()

        prefix = f'{sheet_key}_{cond_label}'

        up_file   = os.path.join(OUT_DIR, f'{prefix}_Up.txt')
        down_file = os.path.join(OUT_DIR, f'{prefix}_Down.txt')

        with open(up_file,   'w') as f: f.write('\n'.join(genes_up))
        with open(down_file, 'w') as f: f.write('\n'.join(genes_down))

        summary.append({
            'Comparison': f'{sheet_key} | {cond_label}',
            'Up (n)':   len(genes_up),
            'Down (n)': len(genes_down),
            'Up file':   os.path.basename(up_file),
            'Down file': os.path.basename(down_file),
        })
        print(f'{prefix}: ↑{len(genes_up)}  ↓{len(genes_down)}')

# Write a summary index
summary_df = pd.DataFrame(summary)
summary_df.to_csv(os.path.join(OUT_DIR, '_summary.csv'), index=False)
print(f'\nAll files saved to: {OUT_DIR}/')
print(f'Total files: {len(summary) * 2} gene lists + 1 summary CSV')
