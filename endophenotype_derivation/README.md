# Endophenotype Derivation — PCA

Derives a PCA-based endophenotype from baseline clinical variables for any cohort. A single script
(`derive_pca.py`) is driven by a per-cohort YAML config. Outputs are named by `cohort_name`.

## Usage

Run from the repository root:

```bash
# CAMP
python endophenotype_derivation/scripts/derive_pca.py \
  --config endophenotype_derivation/config/camp_pca_variables.yaml

# CRA / GACRS
python endophenotype_derivation/scripts/derive_pca.py \
  --config endophenotype_derivation/config/cra_pca_variables.yaml
```

## Dependencies

- Python 3.8+
- `pandas`, `numpy`, `scikit-learn`, `pyyaml`

```bash
pip install pandas numpy scikit-learn pyyaml
```

## YAML config fields

| Field | Required | Description |
|---|---|---|
| `cohort_name` | yes | Used to name all output files |
| `input_file` | yes | Path to phenotype file, resolved relative to repo root |
| `separator` | no | Column delimiter — `"\t"` (default) for TSV, `","` for CSV |
| `output_dir` | no | Output directory; defaults to `endophenotype_derivation/outputs/` |
| `id_column` | no | Subject ID column; falls back to `IID` if present, then row index |
| `variables` | yes | Map of column name → transformation rules (see below) |

## Variable rules

Each key in `variables` is the column name in the input file. Supported fields:

| Field | Required | Description |
|---|---|---|
| `type` | yes | `binary`, `continuous`, or `binarize` |
| `transform` | no | `fraction_to_percent` — multiply by 100 before z-scoring |
| `white_values` | for binarize | listed values → 1, all other non-missing → 0 |
| `low_values` / `high_values` | for binarize | mapped to 0 / 1; stops if any observed value is uncovered |
| `rule: greater_than` + `threshold` | for binarize | `<= threshold` → 0, `> threshold` → 1 |
| `clinical_category` | no | written through to loadings CSV |
| `category_id` | no | written through to loadings CSV |

**binary** — converts to numeric; auto-recodes two unique values to [0, 1] (smaller → 0,
larger → 1) and prints a message; stops with an error if not exactly two unique non-missing values.

**Missing values** — subjects are kept regardless of missingness. After transformation, any missing
value is replaced by the column mean before z-scoring. Per-variable missing counts are printed.

## Outputs

Written to `output_dir` (default `endophenotype_derivation/outputs/`).

### `{cohort_name}_pca_loadings_with_clusters.csv`

One row per feature:

- `variable`, `clinical_category`, `category_id`
- `PC1_loading` … `PC10_loading` (up to 10 PCs)
- `mean_Q1` … `mean_Q5` — mean z-scored value of the feature within each PC1 quintile group
- `max_mean_quintile` — quintile with the highest mean z-score for that feature
- `max_abs_mean_quintile` — quintile with the highest absolute mean z-score

### `{cohort_name}_pca_explained_variance.csv`

One row per PC:

- `PC`, `explained_variance`, `explained_variance_ratio`, `cumulative_explained_variance_ratio`

### `{cohort_name}_pca_subject_scores.csv`

One row per subject (all subjects, since missing values are imputed):

- Subject ID column (name from `id_column`, `IID`, or `subject_id`)
- `PC1_score` … `PC10_score`
- `Cluster` — Q1–Q5 from PC1 quintiles (`pandas.qcut`)
- `numeric_cluster` — 1–5 corresponding to Q1–Q5
