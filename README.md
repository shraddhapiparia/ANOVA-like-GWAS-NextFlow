# ANOVA-like GWAS Pipeline using Nextflow

This repository contains a scalable pipeline for performing genome-wide association testing using ANOVA-style models. It is designed to efficiently process millions of SNPs in batch mode using **Nextflow** and **R**, and supports execution on HPC systems such as SLURM clusters.

---

## Repository Structure

| Path                  | Purpose                                           |
|-----------------------|---------------------------------------------------|
| `workflow/`           | All workflow-related scripts and configs         |
| ├── `new_pipeline.nf` | Main Nextflow workflow definition                |
| ├── `nextflow.config` | Executor settings (e.g., SLURM resources)        |
| ├── `params-file.yaml`| Input paths and batch parameters (template)      |
| ├── `snp_anova.R`     | R script for ANOVA, trend, and linear testing    |
| ├── `nextflow_cra.sh` | SLURM job script to run the pipeline (template)  |
| └── `test.nf`         | Optional test stub for quick validation          |

---

## Requirements

This pipeline requires:

| Tool      | Version  | Notes                              |
|-----------|----------|------------------------------------|
| Java      | 17+      |                                    |
| Nextflow  | 25.10.5+ |                                    |
| R         | 4.3.2    | with `data.table` and `snpStats`   |

The pipeline is designed for SLURM-based HPC. `nextflow.config` sets the executor to `slurm`; change to `local` for local testing.

---

## Required Inputs

| Input | Description |
|-------|-------------|
| PLINK binary fileset | `.bed`, `.bim`, `.fam` files (prefix specified in `params-file.yaml`) |
| Phenotype file | Delimited text file with one row per sample; must contain columns: `Cluster`, `Age`, `Sex`, `PC1`–`PC10`, `pc1_clinical`; sample IDs must match the PLINK `.fam` file |

No toy dataset is included in this repository. You must supply your own PLINK-format genetic data and a matching phenotype file.

---

## Setup

1. Copy and fill in the two template files:

```bash
cp workflow/params-file.yaml workflow/my-params.yaml
# Edit my-params.yaml with your actual paths

cp workflow/nextflow_cra.sh workflow/my-submit.sh
# Edit my-submit.sh with your cluster-specific values (partition, memory, email, paths)
```

2. (Optional) Set a custom R library path in `snp_anova.R` if your HPC requires it.

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
