import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# ── Load & split ──────────────────────────────────────────────────────
df = pd.read_csv('data/csc21_simbad_enriched.csv')
agn     = df[df['source_class'] == 'AGN']
stellar = df[df['source_class'] == 'stellar']

# ── Mann-Whitney U tests ─────────────────────────────────────────────
results = {}
for col in ['hard_hs', 'hard_hm']:
    a = agn[col].dropna()
    s = stellar[col].dropna()
    U, p = mannwhitneyu(a, s, alternative='greater')
    r = 1 - 2 * U / (len(a) * len(s))
    results[col] = dict(U=U, p=p, r=r, n_agn=len(a), n_stel=len(s))
    print(f"{col}:  U = {U:.1f},  p = {p:.4e},  rank-biserial r = {r:.4f}")

# ── Histogram figure ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
bins = np.linspace(-1, 1, 41)

titles = {
    'hard_hs': r'$HR_{HS} = (H-S)/(H+S)$',
    'hard_hm': r'$HR_{HM} = (H-M)/(H+M)$',
}

for ax, col in zip(axes, ['hard_hs', 'hard_hm']):
    a = agn[col].dropna()
    s = stellar[col].dropna()
    r = results[col]

    ax.hist(s, bins=bins, density=True, alpha=0.55, color='#4A90D9',
            label=f'Stellar  (n={r["n_stel"]})', edgecolor='white', linewidth=0.4)
    ax.hist(a, bins=bins, density=True, alpha=0.55, color='#E74C3C',
            label=f'AGN  (n={r["n_agn"]})', edgecolor='white', linewidth=0.4)

    # median lines
    ax.axvline(s.median(), color='#4A90D9', ls='--', lw=1.5, label=f'Stellar median = {s.median():.3f}')
    ax.axvline(a.median(), color='#E74C3C', ls='--', lw=1.5, label=f'AGN median = {a.median():.3f}')

    # annotation
    p_str = f"{r['p']:.2e}" if r['p'] > 1e-30 else f"{r['p']:.1e}"
    ax.text(0.97, 0.95,
            f"Mann-Whitney (AGN > Stellar)\n$U = {r['U']:.0f}$\n$p = {p_str}$\n$r = {r['r']:.3f}$",
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.4', fc='white', alpha=0.85, ec='gray'))

    ax.set_xlabel(titles[col], fontsize=12)
    ax.set_ylabel('Density', fontsize=11)
    ax.legend(fontsize=8.5, loc='upper left')
    ax.set_xlim(-1.05, 1.05)

fig.suptitle('Hardness Ratio Distributions: AGN vs Stellar Sources (CSC 2.1)',
             fontsize=13, fontweight='bold', y=1.02)
fig.tight_layout()
plt.show()
