#DISCLAIMER: This code has been generated using AI, for the prompts used, see: https://claude.ai/share/eac1efb0-d0dc-40cb-a741-7633d6635e2a
"""
Binary classifier: AGN vs. stellar X-ray sources from CSC2.1 + SIMBAD.

Trains on rows labeled 'AGN' or 'stellar', then predicts class probabilities
for (1) unlabeled rows (source_class is NaN) and (2) rows labeled 'other'.
Features: three hardness ratios, log broad-band flux, and |galactic latitude|.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score

INPUT_CSV = "data/csc21_simbad_enriched.csv"
OUTPUT_UNLABELED_CSV = "outputs/20k_unlabeled_predictions.csv"
OUTPUT_OTHER_CSV = "outputs/20k_other_predictions.csv"

FEATURES = ["hard_hs", "hard_hm", "hard_ms", "log_flux_b", "abs_glat"]
RANDOM_STATE = 42


def build_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add engineered columns to df. Does not modify the original."""
    df = df.copy()
    # X-ray fluxes span ~6 orders of magnitude -> log-transform.
    flux = df["flux_aper_b"].astype(float)
    df["log_flux_b"] = np.log10(flux.where(flux > 0, np.nan))
    # Stars cluster toward the Galactic plane; the relevant quantity is |b|.
    df["abs_glat"] = df["glat"].abs()
    return df


def predict_subset(pipe, df_subset: pd.DataFrame, label: str, out_path: str):
    """Run the trained pipeline on a subset, save predictions, print summary."""
    X = df_subset[FEATURES].to_numpy()
    proba = pipe.predict_proba(X)[:, 1]
    pred = np.where(proba >= 0.5, "AGN", "stellar")

    out = df_subset.copy()
    out["prob_AGN"] = proba
    out["predicted_class"] = pred

    out_cols = ["name", "ra", "dec", "glon", "glat",
                "hard_hs", "hard_hm", "hard_ms", "flux_aper_b",
                "simbad_otype", "prob_AGN", "predicted_class"]
    out[out_cols].to_csv(out_path, index=False)

    n_agn = int((pred == "AGN").sum())
    n_stl = int((pred == "stellar").sum())
    print(f"=== Predictions on {label} (n={len(df_subset)}) ===")
    print(f"  AGN:     {n_agn}  ({n_agn / len(df_subset):.1%})")
    print(f"  stellar: {n_stl}  ({n_stl / len(df_subset):.1%})")
    print(f"Wrote: {out_path}\n")


def main():
    df = pd.read_csv(INPUT_CSV)
    df = build_features(df)

    # Three disjoint subsets:
    #   - training: AGN or stellar (have ground-truth labels)
    #   - unlabeled: source_class is NaN (no SIMBAD match at all)
    #   - other: source_class == 'other' (SIMBAD match, but not AGN/stellar)
    train_mask = df["source_class"].isin(["AGN", "stellar"])
    unlabeled_mask = df["source_class"].isna()
    other_mask = df["source_class"] == "other"

    df_train = df.loc[train_mask].copy()
    df_unlabeled = df.loc[unlabeled_mask].copy()
    df_other = df.loc[other_mask].copy()

    print(f"Labeled (AGN+stellar): {len(df_train)}")
    print(f"Unlabeled (NaN):       {len(df_unlabeled)}")
    print(f"'other':               {len(df_other)}")
    print(df_train["source_class"].value_counts().to_string(), "\n")

    X = df_train[FEATURES].to_numpy()
    # Encode AGN=1, stellar=0 so probabilities read as P(AGN).
    y = (df_train["source_class"] == "AGN").astype(int).to_numpy()

    X_tr, X_te, y_tr, y_te = train_test_split(
        X, y, test_size=0.25, stratify=y, random_state=RANDOM_STATE
    )

    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("rf", RandomForestClassifier(
            n_estimators=400,
            max_depth=None,
            min_samples_leaf=2,
            class_weight="balanced",
            n_jobs=-1,
            random_state=RANDOM_STATE,
        )),
    ])
    pipe.fit(X_tr, y_tr)

    # ---- Evaluation on held-out labeled data ----
    y_pred = pipe.predict(X_te)
    y_proba = pipe.predict_proba(X_te)[:, 1]
    print("=== Held-out test performance (1=AGN, 0=stellar) ===")
    print(classification_report(y_te, y_pred, target_names=["stellar", "AGN"], digits=3))
    print("Confusion matrix [rows=true, cols=pred]:")
    print(confusion_matrix(y_te, y_pred))
    print(f"ROC AUC: {roc_auc_score(y_te, y_proba):.4f}\n")

    rf = pipe.named_steps["rf"]
    print("Feature importances (MDI):")
    for name, imp in sorted(zip(FEATURES, rf.feature_importances_),
                            key=lambda t: -t[1]):
        print(f"  {name:12s} {imp:.3f}")
    print()

    # ---- Predict on unlabeled FIRST, then on 'other' ----
    predict_subset(pipe, df_unlabeled, "unlabeled (NaN)", OUTPUT_UNLABELED_CSV)
    predict_subset(pipe, df_other, "'other'", OUTPUT_OTHER_CSV)


if __name__ == "__main__":
    main()
