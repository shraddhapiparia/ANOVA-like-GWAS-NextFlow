#!/usr/bin/env python3
"""
Derive cohort PCA endophenotype from baseline clinical variables.

Reads variable rules from a YAML config, transforms and z-scores features,
runs PCA, derives PC1 quintile group summaries, and writes three CSVs:
  - {cohort_name}_pca_loadings_with_clusters.csv  (feature-level loadings + quintile means)
  - {cohort_name}_pca_explained_variance.csv       (per-PC variance decomposition)
  - {cohort_name}_pca_subject_scores.csv           (subject-level PC scores + cluster labels)

Usage:
    python endophenotype_derivation/scripts/derive_pca.py \
        --config endophenotype_derivation/config/camp_pca_variables.yaml

    python endophenotype_derivation/scripts/derive_pca.py \
        --config endophenotype_derivation/config/cra_pca_variables.yaml
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def binarize_series(col: pd.Series, var_name: str, spec: dict) -> pd.Series:
    col = pd.to_numeric(col, errors="coerce")

    if spec.get("rule") == "greater_than":
        if "threshold" not in spec:
            raise ValueError(f"{var_name}: greater_than rule requires 'threshold'")
        threshold = spec["threshold"]
        print(f"{var_name}: binarized as <= {threshold} -> 0, > {threshold} -> 1")
        return col.map(lambda x: np.nan if pd.isna(x) else (0 if x <= threshold else 1))

    if "white_values" in spec:
        white = set(spec["white_values"])
        return col.map(lambda x: 1 if x in white else (0 if pd.notna(x) else np.nan))

    if "low_values" not in spec or "high_values" not in spec:
        raise ValueError(
            f"{var_name}: binarize requires 'white_values', 'low_values'+'high_values',"
            " or rule: greater_than with 'threshold'"
        )

    low = set(spec["low_values"])
    high = set(spec["high_values"])
    covered = low | high
    observed = set(col.dropna().unique())
    uncovered = observed - covered
    if uncovered:
        raise ValueError(
            f"{var_name}: binarize has values not covered by low_values/high_values: "
            f"{sorted(uncovered)}"
        )

    def _map(x):
        if pd.isna(x):
            return np.nan
        if x in low:
            return 0
        return 1  # must be in high after coverage check

    return col.map(_map)


def transform_variable(col: pd.Series, var_name: str, spec: dict) -> pd.Series:
    vtype = spec["type"]

    if vtype == "binary":
        if "value_map" in spec:
            value_map = spec["value_map"]
            non_missing = col.dropna()
            uncovered = set(non_missing.unique()) - set(value_map.keys())
            if uncovered:
                raise ValueError(
                    f"{var_name}: value_map does not cover all observed non-missing values: "
                    f"{sorted(str(v) for v in uncovered)}"
                )
            col = col.map(lambda x: value_map[x] if pd.notna(x) else np.nan)
        col = pd.to_numeric(col, errors="coerce")
        unique_vals = sorted(col.dropna().unique().tolist())
        if len(unique_vals) == 0:
            raise ValueError(f"{var_name}: binary variable has no non-missing values")
        if len(unique_vals) == 1:
            raise ValueError(
                f"{var_name}: binary variable has only one unique value ({unique_vals[0]})"
                " — constant column cannot be used in PCA"
            )
        if len(unique_vals) > 2:
            raise ValueError(
                f"{var_name}: binary variable has more than two unique values: {unique_vals}"
            )
        if unique_vals != [0, 1]:
            lo, hi = unique_vals
            col = col.map(lambda x: 0 if x == lo else (1 if x == hi else np.nan))
            print(f"{var_name}: binary recoded from {unique_vals} to [0, 1]")

    elif vtype == "continuous":
        col = pd.to_numeric(col, errors="coerce")

    elif vtype == "binarize":
        col = binarize_series(col, var_name, spec)

    else:
        raise ValueError(f"{var_name}: unknown type '{vtype}'")

    transform = spec.get("transform")
    if transform == "fraction_to_percent":
        col = col * 100
    elif transform is not None:
        raise ValueError(f"{var_name}: unknown transform '{transform}'")

    return col


def main():
    parser = argparse.ArgumentParser(
        description="Derive cohort PCA endophenotype from baseline clinical variables."
    )
    parser.add_argument("--config", required=True, type=Path, help="Path to YAML config file")
    args = parser.parse_args()

    config_path = args.config.resolve()
    if not config_path.exists():
        sys.exit(f"ERROR: config file not found: {config_path}")

    config = load_config(config_path)

    # input_file resolved relative to CWD (expected: repo root)
    input_file = Path(config["input_file"])
    if not input_file.is_absolute():
        input_file = Path.cwd() / input_file
    if not input_file.exists():
        sys.exit(f"ERROR: input file not found: {input_file}")

    output_dir = Path(config.get("output_dir", config_path.parent.parent / "outputs"))
    if not output_dir.is_absolute():
        output_dir = Path.cwd() / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_name = config.get("cohort_name", "cohort")

    variables = config.get("variables")
    if not variables:
        sys.exit("ERROR: 'variables' section is empty or missing in config")

    # Load input — separator defaults to tab; set separator in YAML for CSV inputs
    sep = config.get("separator", "\t")
    df = pd.read_csv(input_file, sep=sep, low_memory=False)
    print(f"Loaded {len(df)} subjects from {input_file.name}")

    # Transform each variable
    transformed = {}
    for var_name, spec in variables.items():
        if var_name not in df.columns:
            sys.exit(f"ERROR: column '{var_name}' not found in input file")
        transformed[var_name] = transform_variable(df[var_name].copy(), var_name, spec)

    pca_df = pd.DataFrame(transformed, index=df.index)

    # Check for entirely missing variables before imputation
    all_missing = pca_df.columns[pca_df.isna().all()].tolist()
    if all_missing:
        sys.exit(f"ERROR: the following variables are entirely missing and cannot be imputed: {all_missing}")

    # Report missing counts and mean-impute
    n_subjects = len(pca_df)
    n_any_missing = int(pca_df.isna().any(axis=1).sum())
    missing_counts = pca_df.isna().sum()
    missing_counts = missing_counts[missing_counts > 0]

    print(f"Subjects total: {n_subjects}")
    print(f"Subjects with at least one missing PCA variable: {n_any_missing}")
    if missing_counts.empty:
        print("No missing values found.")
    else:
        print("Per-variable missing counts:")
        for var_name, count in missing_counts.items():
            print(f"  {var_name}: {count} missing")
        col_means = pca_df.mean()
        for col in missing_counts.index:
            pca_df[col] = pca_df[col].fillna(col_means[col])
        print("Missing values mean-imputed before z-scoring.")

    feature_names = list(pca_df.columns)

    # Z-score all features
    scaler = StandardScaler()
    X = scaler.fit_transform(pca_df.values)  # shape: (n_subjects, n_features)

    # PCA (up to 10 components)
    n_components = min(10, len(feature_names), len(pca_df))
    pca = PCA(n_components=n_components, random_state=42)
    scores = pca.fit_transform(X)  # shape: (n_subjects, n_components)

    # PC1 quintiles
    pc1 = pd.Series(scores[:, 0])
    quintiles = pd.qcut(pc1, q=5, labels=["Q1", "Q2", "Q3", "Q4", "Q5"])

    # Mean z-scored value per feature per quintile group
    X_df = pd.DataFrame(X, columns=feature_names)
    X_df["_quintile"] = quintiles.values
    quintile_means = {
        f"mean_{q}": X_df.loc[X_df["_quintile"] == q, feature_names].mean()
        for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]
    }

    # Build feature-level loadings output
    rows = []
    for i, var_name in enumerate(feature_names):
        spec = variables[var_name]
        row = {
            "variable": var_name,
            "clinical_category": spec.get("clinical_category"),
            "category_id": spec.get("category_id"),
        }
        for pc_idx in range(n_components):
            row[f"PC{pc_idx + 1}_loading"] = pca.components_[pc_idx, i]
        for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
            row[f"mean_{q}"] = quintile_means[f"mean_{q}"][var_name]
        rows.append(row)

    loadings_df = pd.DataFrame(rows)

    mean_cols = [f"mean_{q}" for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]]
    means = loadings_df[mean_cols]
    loadings_df["max_mean_quintile"] = means.idxmax(axis=1).str.replace("mean_", "", regex=False)
    loadings_df["max_abs_mean_quintile"] = means.abs().idxmax(axis=1).str.replace("mean_", "", regex=False)

    loadings_out = output_dir / f"{cohort_name}_pca_loadings_with_clusters.csv"
    loadings_df.to_csv(loadings_out, index=False)
    print(f"Saved: {loadings_out}")

    # Variance decomposition output
    variance_df = pd.DataFrame({
        "PC": [f"PC{i + 1}" for i in range(n_components)],
        "explained_variance": pca.explained_variance_,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "cumulative_explained_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
    })

    variance_out = output_dir / f"{cohort_name}_pca_explained_variance.csv"
    variance_df.to_csv(variance_out, index=False)
    print(f"Saved: {variance_out}")

    # Subject-level output: resolve ID column
    id_col_name = config.get("id_column")
    if id_col_name and id_col_name in df.columns:
        id_label = id_col_name
        id_values = df[id_col_name].values
    elif "IID" in df.columns:
        id_label = "IID"
        id_values = df["IID"].values
    else:
        id_label = "subject_id"
        id_values = df.index.values

    subject_data = {id_label: id_values}
    for pc_idx in range(n_components):
        subject_data[f"PC{pc_idx + 1}_score"] = scores[:, pc_idx]
    subject_data["Cluster"] = quintiles.values
    subject_data["numeric_cluster"] = (
        pd.Series(quintiles).map({"Q1": 1, "Q2": 2, "Q3": 3, "Q4": 4, "Q5": 5}).values
    )

    subject_df = pd.DataFrame(subject_data)
    subject_out = output_dir / f"{cohort_name}_pca_subject_scores.csv"
    subject_df.to_csv(subject_out, index=False)
    print(f"Saved: {subject_out}")


if __name__ == "__main__":
    main()
