# ANOVA-like GWAS Pipeline using Nextflow

This repository contains a scalable pipeline for performing genome-wide association testing using ANOVA-style models. It is designed to efficiently process millions of SNPs in batch mode using **Nextflow** and **R**, and supports execution on HPC systems such as SLURM clusters.

---

## 📁 Repository Structure

| Path                  | Purpose                                           |
|-----------------------|---------------------------------------------------|
| `workflow/`           | All workflow-related scripts and configs         |
| ├── `new_pipeline.nf` | Main Nextflow workflow definition                |
| ├── `nextflow.config` | Executor settings (e.g., SLURM resources)        |
| ├── `params-file.yaml`| Input paths and batch parameters                 |
| ├── `snp_anova.R`     | R script for ANOVA, trend, and linear testing    |
| ├── `nextflow_cra.sh` | SLURM job script to run the pipeline             |
| └── `test.nf`         | Optional test stub for quick validation          |
| `data/`               | Input PLINK data folder (e.g., `.bed/.bim/.fam`) |
| `results/`            | Output directory (auto-populated)                |
| `docs/`               | Project notes, figures, or plots                 |

---

## ⚙️ Requirements

This pipeline is built using **Nextflow DSL2** and requires the following environment setup:

| Tool       | Version      |
|------------|--------------|
| Java       | 17+           |
| R          | 4.3.2        |
| Nextflow   | 25.10.5      |

---

## ▶️ How to Run

### On SLURM-based HPC
Submit with:

cd workflow

sbatch nextflow_cra.sh

Or run directly

cd workflow

nextflow run new_pipeline.nf -params-file params-file.yaml -resume

## 📊 Output

Each batch produces a result file containing:

- SNP metadata (ID, chromosome, position)
- ANOVA F-statistic and p-value
- Trend-based association coefficients
- Correlation with top principal component

Final results can be aggregated for:

- Manhattan plots
- Q-Q plots
- Endotype-stratified association signals

