import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/csc21_simbad_enriched.csv')

stellar = df[df['source_class'] == 'stellar']
agn = df[df['source_class'] == 'AGN']

fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(stellar['ra'], stellar['dec'], s=10, alpha=0.6, color='tomato', label='Stellar', marker='*')
ax.scatter(agn['ra'], agn['dec'], s=10, alpha=0.6, color='steelblue', label='AGN', marker='D')

ax.set_xlabel('Right Ascension (deg)')
ax.set_ylabel('Declination (deg)')
ax.set_title('CSC 2.1 Sources: Stellar vs AGN')
ax.legend()
ax.invert_xaxis()  # RA increases to the left by convention

plt.tight_layout()
plt.savefig('scatter_stellar_agn.png', dpi=150)
plt.show()
