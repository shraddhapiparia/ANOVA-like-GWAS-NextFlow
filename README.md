# ANOVA-like GWAS Pipeline using Nextflow

This repository contains a scalable pipeline for performing genome-wide association testing using ANOVA-style models. It is designed to efficiently process millions of SNPs in batch mode using **Nextflow** and **R**, and supports execution on HPC systems such as SLURM clusters.

---

## Repository Structure

| Path                  | Purpose                                           |
|-----------------------|---------------------------------------------------|
| `workflow/`            | All workflow-related scripts and configs          |
| ├── `new_pipeline.nf`  | Main Nextflow workflow definition                 |
| ├── `test.nf`          | Smoke-test workflow (no module load, uses example/) |
| ├── `nextflow.config`  | Executor settings (e.g., SLURM resources)         |
| ├── `params-file.yaml` | Input paths and batch parameters (template)       |
| ├── `params-test.yaml` | Params for smoke test (fill in paths)             |
| ├── `snp_anova.R`      | R script for ANOVA, trend, and linear testing     |
| └── `nextflow_cra.sh`  | SLURM job script to run the pipeline (template)   |
| `example/`             | Synthetic dataset for smoke-testing               |
| ├── `example.bed/bim/fam` | Tiny PLINK fileset (25 samples, 10 SNPs)       |
| ├── `example_pheno.csv`| Matching phenotype file                           |
| └── `make_example_data.py` | Script that regenerates the example files    |

---

## Dependencies

This pipeline was developed and run on an HPC cluster where software is provided via `module load`. The pipeline assumes the following are available:

| Software | Version used | Notes |
|----------|-------------|-------|
| Java | 17+ | Required by Nextflow |
| Nextflow | 25.10.5+ | DSL2 |
| R | 4.3.2 | |
| R: `data.table` | any recent | CRAN |
| R: `snpStats` | any recent | Bioconductor |

See `workflow/nextflow_cra.sh` for the exact `module load` commands used in production. There is no conda environment or container — dependencies are assumed to be pre-installed on the cluster.

**Smoke test (local):** If running `test.nf` locally rather than on HPC, you need `Rscript` on your PATH with `data.table` and `snpStats` installed, plus Java and Nextflow. Install `snpStats` via:
```r
install.packages("BiocManager"); BiocManager::install("snpStats")
```

---

## Required Inputs

| Input | Description |
|-------|-------------|
| PLINK binary fileset | `.bed`, `.bim`, `.fam` files (prefix specified in `params-file.yaml`) |
| Phenotype file | Delimited text file with one row per sample; must contain columns: `Cluster`, `Age`, `Sex`, `PC1`–`PC10`, `pc1_clinical`; sample IDs must match the PLINK `.fam` file |

A tiny synthetic dataset (25 samples, 10 SNPs) is included in `example/` for smoke-testing. It is not real data and exists only to verify the pipeline runs end-to-end. For real analyses you must supply your own PLINK-format genetic data and a matching phenotype file.

---

## Smoke Test (example data)

Verifies that Nextflow, R, `data.table`, and `snpStats` are all wired up correctly.
Runs `snp_anova.R` on a synthetic 25-sample / 10-SNP dataset without any HPC modules.

**Requires:** Java 17+, Nextflow, R ≥ 4.3.2 with `data.table` and `snpStats`, `Rscript` on `PATH`.
**Does not require:** SLURM, `module load`, or real genetic data.

```bash
# 1. Edit params-test.yaml — replace /absolute/path/to/repo with your actual path
#    (Three fields: plink_prefix, pheno_file, r_script)

# 2. Run from the workflow/ directory
cd workflow
nextflow run test.nf -params-file params-test.yaml
```

Expected output in `workflow/` (or the Nextflow work directory):

```
gwas_results_batch1.csv   # 10 rows, one per SNP
```

The CSV will contain columns `snp`, `chromosome`, `position`, `ANOVA_F`, `ANOVA_pvalue`,
`Trend_Estimate`, `Trend_pvalue`, `PC1_Clinical_Estimate`, `PC1_Clinical_pvalue`.
Statistical values are not meaningful (synthetic data); the test passes if the file is produced with 10 rows and no R errors.

---

## Pipeline Setup

Once the environment is ready (see above), configure the pipeline:

```bash
cp workflow/params-file.yaml workflow/my-params.yaml
# Edit my-params.yaml — fill in plink_prefix, pheno_file, r_script, bim_file

cp workflow/nextflow_cra.sh workflow/my-submit.sh
# Edit my-submit.sh — fill in partition, memory, email, log paths, TMPDIR
```

If your HPC requires a non-default R library path, uncomment and set `.libPaths(...)` at line 7 of `workflow/snp_anova.R`.

---

## How to Run

### On SLURM-based HPC

```bash
cd workflow
sbatch my-submit.sh
```

### Directly (interactive or local testing)

```bash
cd workflow
nextflow run new_pipeline.nf -params-file my-params.yaml -resume
```

---

## Expected Outputs

Each batch produces one CSV file: `gwas_results_batch{N}.csv`

Each file contains one row per SNP with columns:

| Column | Description |
|--------|-------------|
| `snp` | SNP identifier |
| `chromosome`, `position`, `allele.1`, `allele.2` | SNP metadata |
| `ANOVA_F`, `ANOVA_pvalue` | ANOVA F-statistic and p-value (cluster as factor) |
| `Trend_Estimate`, `Trend_pvalue` | Linear trend test (ordinal cluster encoding) |
| `PC1_Clinical_Estimate`, `PC1_Clinical_pvalue` | Linear regression vs. clinical PC1 |

Batch outputs can be concatenated and used to generate Manhattan plots, Q-Q plots, and endotype-stratified association signals.

---

## Cluster-specific vs. Portable

**Cluster-specific** (must be adapted per site):
- `nextflow_cra.sh` — SLURM directives, module names, paths
- `nextflow.config` — executor, memory/time per process
- `snp_anova.R` line 7 — custom R library path (`.libPaths`)
- `new_pipeline.nf` — `module load R/4.3.2` and `TMPDIR` inside the `RUN_BATCH` script block

**Portable** (no changes needed):
- `new_pipeline.nf` — workflow logic
- `snp_anova.R` — statistical analysis
- `params-file.yaml` — parameter structure (only paths need filling in)
