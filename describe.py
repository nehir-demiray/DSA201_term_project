import pandas as pd

pd.set_option('display.float_format', '{:.4e}'.format)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 120)

df = pd.read_csv('data/csc21_simbad_enriched.csv')

cols = ['hard_hs', 'hard_hm', 'hard_ms', 'flux_aper_b', 'flux_aper_s', 'flux_aper_m', 'flux_aper_h']

df['source_class'] = df['source_class'].fillna('unmatched')

for cls in ['stellar', 'AGN', 'other', 'unmatched']:
    subset = df[df['source_class'] == cls]
    print(f'\n=== {cls} (n={len(subset)}) ===')
    print(subset[cols].describe().to_string())
