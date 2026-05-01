# Endotype-Aware GWAS via Categorical ANCOVA
This repository contains a scalable pipeline for genome-wide association testing across PCA-derived asthma endophenotype groups. It includes:

1. A Python module to derive clinical PCA endophenotype quintiles from cohort phenotype files.
2. A Nextflow/R workflow to test SNP associations across those endophenotype groups using ANOVA-style models.

The workflow is designed for batch processing millions of SNPs on HPC systems such as SLURM clusters.

---

## Repository Structure

| Path | Purpose |
|------|---------|
| `endophenotype_derivation/` | Derives PCA-based clinical endophenotype quintiles |
| `endophenotype_derivation/config/` | YAML configs for cohort-specific PCA variables |
| `endophenotype_derivation/scripts/derive_pca.py` | Cohort-agnostic PCA/endophenotype derivation script |
| `workflow/` | Nextflow/R GWAS workflow |
| `workflow/main.nf` | Main Nextflow workflow |
| `workflow/test.nf` | Smoke-test workflow using synthetic example data |
| `workflow/snp_anova.R` | SNP-level ANOVA, trend, and clinical PC1 association testing |
| `workflow/nextflow.config` | SLURM profile and resource settings |
| `workflow/params-file.yaml` | Input paths and batch parameters |
| `workflow/nextflow_cra.sh` | SLURM submission script template |
| `example/` | Synthetic PLINK dataset for smoke testing |

---

## Endophenotype Derivation

The `endophenotype_derivation/` module derives PCA-based clinical endophenotype groups from cohort phenotype files.

For each cohort, selected baseline clinical variables are:

1. transformed according to a YAML config,
2. mean-imputed if missing,
3. z-scored,
4. summarized using PCA,
5. assigned to PC1 quintiles Q1–Q5.

These quintiles represent ordered asthma endophenotype groups along the dominant clinical variation axis.

### Run CAMP PCA

```bash
python endophenotype_derivation/scripts/derive_pca.py \
  --config endophenotype_derivation/config/camp_pca_variables.yaml
```

### Run CRA/GACRS PCA

```bash
python endophenotype_derivation/scripts/derive_pca.py \
  --config endophenotype_derivation/config/cra_pca_variables.yaml
```

## PCA Outputs

Outputs are written to `endophenotype_derivation/outputs/` and are not tracked by git.

| Output | Description |
|--------|-------------|
| `{cohort}_pca_loadings_with_clusters.csv` | Feature-level PCA loadings plus Q1–Q5 feature summaries |
| `{cohort}_pca_explained_variance.csv` | Per-PC explained variance |
| `{cohort}_pca_subject_scores.csv` | Subject-level PC scores and Q1–Q5 endophenotype labels |

The subject score file can be used to create the phenotype input for the downstream ANOVA GWAS workflow.

---

## Key Results

## Key Results

This pipeline produced the methodology and findings reported in two publications: the **PCA-based endophenotype derivation** (*Respiratory Research* 2025) and the **ANCOVA-based GWAS** (*J. Pers. Med.* 2026).

### Endophenotype Derivation (PCA)

PCA was applied to ~20 standardized baseline clinical features (demographics, lung function, atopic status, exacerbation history) across three independent pediatric asthma cohorts: CAMP (n=1,041), PACT (n=230), and GACRS (n=1,165). Subjects were assigned to ordinal endophenotypes Q1–Q5 based on PC1 score quintiles.

- **Reproducibility across cohorts:** Atopy, lung function, and demographics were the dominant contributors to PC1 in all three cohorts (CAMP: 67%, PACT: 49%, GACRS: 60% of total loading magnitude).
- **Severity gradient:** Pre-bronchodilator FEV1% predicted declined monotonically from Q1 → Q5, while IgE, eosinophil counts, and SABA usage increased — consistent across cohorts.
- **Treatment-response prediction:** In CAMP, ICS-treated participants showed a progressive 1-year FEV1 gain from Q1 (median 0.0%) to Q5 (median 16.7%), versus a much smaller gradient on placebo (–2.4% → 5.5%). PACT replicated this pattern with fluticasone vs montelukast.

### Endotype-Aware GWAS (ANCOVA)

For each SNP, the pipeline fits:
``` genotype ~ endophenotype + age + sex + PC1–PC10```

where `endophenotype` is a 5-level categorical factor (Q1–Q5) capturing severity gradient. Three statistics are reported per SNP:

- **ANOVA F-test**: omnibus test for any allele-frequency difference across endophenotypes
- **Trend test**: linear allele-frequency gradient across ordinal endophenotypes
- **PC1-clinical regression**: continuous severity association

Significant SNPs from the ANOVA F-test are followed up with Tukey HSD post-hoc contrasts to identify which subtype pairs drive the signal.

**Findings in CAMP (discovery, n=792) and GACRS (replication, n=1,030):**

- 244 genome-wide significant SNPs in CAMP (p ≤ 5×10⁻⁸); 6 LD-independent loci after clumping
- All 6 loci nominally replicated in GACRS — **rs10964536, rs28892326, rs2823880, rs10086065, rs12448208, rs2754324**
- Risk-allele frequency increased monotonically with severity in 5 of 6 loci across both cohorts (delta MAF 7.5–12% in CAMP, ≥4% in GACRS)
- ANCOVA recovered **13 significant post-hoc contrasts in CAMP** vs only **4 surviving Bonferroni-corrected one-vs-rest logistic regression** — confirming improved power over the conventional pairwise approach
- Several loci map to genes with airway-relevant biology (e.g., *DGKI* — airway smooth muscle remodeling; *MIR99AHG* — host gene of the miR-99a/let-7c/miR-125b-2 cluster).

---

## Where This Fits in the Genomics Workflow

This pipeline operates on **PLINK-format genotypes downstream of standard NGS preprocessing** (QC, alignment, variant calling, imputation). Specifically, it expects:

- Imputed, QC'd PLINK binary filesets (`.bed`/`.bim`/`.fam`)
- A phenotype file with PCA-derived endophenotype labels and ancestry covariates

For an **end-to-end example covering FASTQ → variants**, see the companion repo:
[from-fastq-to-asthma-gwas](https://github.com/shraddhapiparia/from-fastq-to-asthma-gwas)

For the **PCA-based endotype derivation** step that produces the `Cluster` and `pc1_clinical` columns expected by this pipeline, see [`endotype_derivation/`](./endotype_derivation/) (or the *Respiratory Research* 2025 paper, Supplementary Table S1, for the loading matrix to apply to your own clinical data).

---

## Dependencies

This pipeline was developed for an HPC environment where software is provided via `module load`.

| Software | Version used | Notes |
|----------|--------------|-------|
| Java | 17+ | Required by Nextflow |
| Nextflow | 25.10.5+ | DSL2 |
| R | 4.3.2 | Used for SNP-level association testing |
| R: `data.table` | recent | CRAN |
| R: `snpStats` | recent | Bioconductor |
| Python | 3.10+ | Used for PCA endophenotype derivation |
| Python: `pandas`, `numpy`, `scikit-learn`, `pyyaml` | recent | PCA module |

Install `snpStats` in R with:

```r
install.packages("BiocManager")
BiocManager::install("snpStats")
```

## Required Inputs for GWAS Workflow

| Input | Description |
|-------|-------------|
| PLINK binary fileset | .bed, .bim, .fam files |
| Phenotype file | Delimited file with one row per sample |

The phenotype file must contain:

- Cluster
- Age
- Sex
- PC1-PC10
- pc1_clinical

Cluster should represent the PCA-derived Q1-Q5 endophenotype label.  
pc1_clinical should represent the continuous clinical PC1 score.

The phenotype file must also contain sample IDs matching the PLINK .fam file. The R script aligns phenotype rows to PLINK sample order by sample ID and stops if required samples are missing.

---

## Smoke Test

A tiny synthetic dataset is included in example/ to verify that Nextflow, R, data.table, and snpStats are wired correctly.

    cd workflow
    nextflow run test.nf -params-file params-test.yaml

Expected output:

    gwas_results_batch1.csv

The file should contain 10 rows, one per synthetic SNP. Statistical values are not meaningful.

---

## Pipeline Setup

Copy and edit the parameter template:

    cp workflow/params-file.yaml workflow/my-params.yaml

Fill in:

- plink_prefix
- pheno_file
- r_script
- bim_file

If your HPC requires a custom R library path, update .libPaths(...) in workflow/snp_anova.R.

---

## How to Run

### SLURM

    cd workflow
    sbatch nextflow_cra.sh

### Interactive HPC session

    cd workflow
    nextflow run main.nf -params-file my-params.yaml -profile slurm -resume

---

## GWAS Outputs

The workflow partitions SNPs into independent batches. Each batch writes:

    gwas_results_batch{N}.csv

Each row represents one SNP.

| Column | Description |
|--------|-------------|
| snp | SNP identifier |
| chromosome, position, allele.1, allele.2 | SNP metadata |
| ANOVA_F, ANOVA_pvalue | Association with endophenotype group as a factor |
| Trend_Estimate, Trend_pvalue | Ordinal Q1-Q5 trend test |
| PC1_Clinical_Estimate, PC1_Clinical_pvalue | Association with continuous clinical PC1 |

Significant ANOVA hits can be followed up separately with Tukey HSD post-hoc contrasts to identify which endophenotype groups differ.


