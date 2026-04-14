#!/usr/bin/env python
"""
Query Chandra Source Catalog (CSC 2.1) for the top N X-ray sources by
detection significance, cross-match with SIMBAD for object-type labels,
and save the enriched catalog to CSV.

Requirements:
    pip install pyvo astroquery astropy numpy pandas
"""

import warnings
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.simbad import Simbad
import pyvo as vo

warnings.filterwarnings("ignore", category=UserWarning)

# ── Configuration ──────────────────────────────────────────────
N_SOURCES        = 5000    # top N by detection significance
MIN_SIGNIFICANCE = 3.0    # minimum detection significance
XMATCH_RADIUS    = 3.0    # SIMBAD cross-match radius (arcsec)
OUTPUT_CSV       = "data/csc21_simbad_enriched.csv"

# ── Classification dictionaries ────────────────────────────────
AGN_TYPES = {
    "AGN", "AG?", "Sy1", "Sy2", "Bla", "BLL", "QSO", "Q?",
    "LIN", "SyG", "rG",  "H2G", "LSB", "SBG", "bCG", "GiG",
    "GiP", "EmG", "BiC", "G",   "GNe", "GiC", "Cl*G",
}

STELLAR_TYPES = {
    "*",   "**",  "HB*", "RGB*", "AGB*", "WD*", "Psr", "XB*",
    "LXB", "HXB", "CV*", "V*",  "PM*",  "Em*", "Be*", "ro*",
    "TT*", "Ce*", "RR*", "Mi*", "sg*",  "s*r", "s*y", "s*b",
    "WR*", "N*",  "gam", "SNR", "Pla",  "RS*", "BY*", "El*",
    "SX*", "No*", "dS*", "RG*", "EB*",  "Al*", "bL*", "Y*O",
    "pr*", "BS*", "OH*", "CH*", "MS*",  "S*",  "LP*", "AB*",
}


def classify(otype):
    if otype is None:
        return None
    otype = str(otype).strip()
    if otype in AGN_TYPES:
        return "AGN"
    if otype in STELLAR_TYPES:
        return "stellar"
    return "other"


# ==============================================================
# 1. Query CSC 2.1
# ==============================================================
print(f"[1/4] Querying CSC 2.1 for top {N_SOURCES} sources ...")

tap = vo.dal.TAPService("http://cda.cfa.harvard.edu/csc21tap")

query = f"""
SELECT TOP {N_SOURCES}
    m.name,
    m.ra,
    m.dec,
    m.err_ellipse_r0,
    m.significance,
    m.hard_hs,
    m.hard_hm,
    m.hard_ms,
    m.flux_aper_b,
    m.flux_aper_s,
    m.flux_aper_m,
    m.flux_aper_h,
    m.extent_flag,
    m.var_flag
FROM csc21.master_source m
WHERE m.significance > {MIN_SIGNIFICANCE}
  AND m.hard_hs IS NOT NULL
ORDER BY m.significance DESC
"""

result = tap.search(query)
csc = result.to_table()
print(f"   Retrieved {len(csc)} sources.")

# Add galactic coordinates
coords = SkyCoord(ra=csc["ra"], dec=csc["dec"], unit="deg", frame="icrs")
gal = coords.galactic
csc["glon"] = gal.l.deg
csc["glat"] = gal.b.deg

# ==============================================================
# 2. Cross-match with SIMBAD
# ==============================================================
print(f"[2/4] Cross-matching with SIMBAD (r = {XMATCH_RADIUS}″) ...")

simbad = Simbad()
simbad.add_votable_fields("otype")
simbad.ROW_LIMIT = 1

otypes, simbad_ids = [], []
for i, coord in enumerate(coords):
    if (i + 1) % 100 == 0 or (i + 1) == len(coords):
        print(f"      ... {i+1}/{len(coords)}")
    try:
        r = simbad.query_region(coord, radius=XMATCH_RADIUS * u.arcsec)
        if r is not None and len(r) > 0:
            otypes.append(r["otype"][0])
            simbad_ids.append(r["main_id"][0])
        else:
            otypes.append(None)
            simbad_ids.append(None)
    except Exception:
        otypes.append(None)
        simbad_ids.append(None)

csc["simbad_otype"] = otypes
csc["simbad_id"]    = simbad_ids

n_matched = sum(1 for o in otypes if o is not None)
print(f"   Matched {n_matched}/{len(csc)} sources.")

# ==============================================================
# 3. Classify
# ==============================================================
print("[3/4] Classifying sources ...")
csc["source_class"] = [classify(o) for o in otypes]

# ==============================================================
# 4. Convert to pandas and save
# ==============================================================
df = csc.to_pandas()

# Clean up masked/object columns to proper numeric types
for col in ["hard_hs", "hard_hm", "hard_ms",
            "flux_aper_b", "flux_aper_s", "flux_aper_m", "flux_aper_h",
            "err_ellipse_r0", "significance"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df.to_csv(OUTPUT_CSV, index=False)
print(f"[4/4] Saved {len(df)} rows to {OUTPUT_CSV}")

# Summary
counts = df["source_class"].value_counts(dropna=False)
print("\n   Class breakdown:")
for cls, n in counts.items():
    label = cls if pd.notna(cls) else "unmatched"
    print(f"      {label:12s}: {n}")
