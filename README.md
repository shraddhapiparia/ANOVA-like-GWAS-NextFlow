# Endotype-Aware GWAS via Categorical ANCOVA

A reproducible Nextflow pipeline for detecting genetic variants associated with clinically defined disease subtypes (endophenotypes), rather than traditional case-control labels. By modeling all subtypes simultaneously using ANCOVA with Tukey post-hoc contrasts, this approach identifies subtype-specific signals that are missed by one-vs-rest logistic regression — published in *J. Pers. Med.* (2026), where it identified 244 genome-wide significant SNPs in CAMP, with six LD-independent loci replicating in GACRS.

**Why ANCOVA over case-control GWAS?** Conventional GWAS treats all cases as one group, implicitly assuming shared genetic architecture across heterogeneous disease subtypes. One-vs-rest logistic regression, the common alternative, multiplies the testing burden and loses power. ANCOVA tests all subtypes within a single model, recovers shared and subtype-specific signal in one pass, and reduces multiple-testing inflation.

---

## Repository Structure

| Path                  | Purpose                                           |
|-----------------------|---------------------------------------------------|
| `workflow/`            | All workflow-related scripts and configs          |
| ├── `main.nf`          | Main Nextflow workflow definition                 |
| ├── `test.nf`          | Smoke-test workflow (no module load, uses example/) |
| ├── `nextflow.config`  | Profiles and resource settings (slurm profile)    |
| ├── `params-file.yaml` | Input paths and batch parameters (template)       |
| ├── `params-test.yaml` | Params for smoke test (fill in paths)             |
| ├── `snp_anova.R`      | R script for ANOVA, trend, and linear testing     |
| └── `nextflow_cra.sh`  | SLURM job script to run the pipeline (template)   |
| `example/`             | Synthetic dataset for smoke-testing               |
| ├── `example.bed/bim/fam` | Tiny PLINK fileset (25 samples, 10 SNPs)       |
| ├── `example_pheno.csv`| Matching phenotype file                           |
| └── `make_example_data.py` | Script that regenerates the example files    |

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

Sample order must match the PLINK .fam file row order; the current implementation does not merge by sample ID.

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

### Directly on HPC (interactive session)

```bash
cd workflow
nextflow run main.nf -params-file my-params.yaml -profile slurm -resume
```

---

## Expected Outputs


The pipeline parallelizes association testing by partitioning genome-wide SNPs into independent batches processed concurrently by Nextflow. Each batch produces one CSV file: `gwas_results_batch{N}.csv`

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
- `nextflow_cra.sh` — SLURM directives, log paths, node/partition names
- `nextflow.config` `slurm` profile — `module load` names, `TMPDIR`, memory/time limits
- `snp_anova.R` line 7 — optional custom R library path (`.libPaths`)

**Portable** (no changes needed):
- `main.nf` — workflow logic and process definitions
- `snp_anova.R` — statistical analysis
- `params-file.yaml` — parameter structure (only paths need filling in)

---

## Publications

This repository implements the methodology described in:

- **Piparia S, Hadikhani P, Ziniti J, Hecker J, Kho A, Sharma R, Celedón JC, Weiss ST, McGeachie M, Tantisira K.** *A Categorical ANCOVA Approach to Severity Endophenotype-Specific GWAS in Childhood Asthma.* J. Pers. Med. 16(1):32 (2026). [DOI](https://www.mdpi.com/2075-4426/16/1/32)
- **Piparia S, Kho A, Desai B, Wong R, Sharma R, Celedón JC, Weiss ST, McGeachie M, Tantisira K.** *A principal component analysis-based endophenotype definition for change in lung function and inhaled corticosteroid treatment response in childhood asthma.* Respiratory Research 26:351 (2025). [DOI](https://link.springer.com/article/10.1186/s12931-025-03426-z)

The 2025 *Respiratory Research* paper defines the PC1-based severity endophenotypes (Q1–Q5) used as input here. The 2026 *J. Pers. Med.* paper applies the ANCOVA framework in this repository to identify subtype-specific genetic variants in CAMP (n=792) and replicates them in GACRS (n=1,030).

