args <- commandArgs(trailingOnly=TRUE)
plink_prefix <- args[1]
pheno_file <- args[2]
start_idx <- as.numeric(args[3])
end_idx <- as.numeric(args[4])
batch_num <- as.numeric(args[5])
library(data.table)
library(snpStats)

# Optional: set custom R library path for HPC environments where
# Bioconductor packages live outside the default location.
# Example:
#   .libPaths(c("/cluster/r-libs/4.3.2", .libPaths()))

phenotype_data <- fread(pheno_file)

# --- Validate phenotype file has all required columns -----------------------
required_cols <- c("IID", "Cluster", "Age", "Sex",
                   paste0("PC", 1:10), "pc1_clinical")
missing_cols <- setdiff(required_cols, names(phenotype_data))
if (length(missing_cols) > 0) {
  stop("Phenotype file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

phenotype_data$IID <- as.character(phenotype_data$IID)
phenotype_data$Cluster <- as.factor(phenotype_data$Cluster)
phenotype_data$numeric_cluster <- as.numeric(phenotype_data$Cluster)

# --- Read SNP metadata from .bim --------------------------------------------
bim_file <- paste0(plink_prefix, ".bim")
snp_info <- fread(bim_file, header = FALSE)
colnames(snp_info) <- c("chromosome", "snp_id", "cM", "position",
                        "allele.1", "allele.2")
snp_info <- snp_info[!duplicated(snp_id)]
snp_ids <- snp_info$snp_id
num_snps <- length(snp_ids)

snps_batch_ids <- snp_ids[start_idx:end_idx]
snp_info_batch <- snp_info[match(snps_batch_ids, snp_info$snp_id), ]

# --- Read PLINK sample order from .fam --------------------------------------
# Read .fam directly rather than relying on snpStats row names: snpStats
# picks pedigree IDs when unique and falls back to member IDs otherwise,
# so reading column 2 ourselves is unambiguous.
fam_file <- paste0(plink_prefix, ".fam")
fam_data <- fread(fam_file, header = FALSE,
                  col.names = c("FID", "IID", "PAT", "MAT", "Sex_fam", "Pheno_fam"))
fam_iids <- as.character(fam_data$IID)

# --- Read genotypes for this batch ------------------------------------------
cat("\nReading batch", batch_num, ":", start_idx, "-", end_idx, "\n")
genetic_data_batch <- read.plink(plink_prefix, select.snps = snps_batch_ids)
genotypes_batch <- as(genetic_data_batch$genotypes, "numeric")

# --- Align phenotype rows to PLINK sample order by IID ----------------------
# We do not assume the phenotype file is pre-sorted to match the .fam file.
# Every .fam IID must be present in the phenotype file; extra phenotype rows
# (samples not in the genotype data) are dropped.
pheno_idx <- match(fam_iids, phenotype_data$IID)
n_missing <- sum(is.na(pheno_idx))
if (n_missing > 0) {
  missing_iids <- fam_iids[is.na(pheno_idx)]
  stop(n_missing, " PLINK sample(s) missing from phenotype file. ",
       "First few: ",
       paste(head(missing_iids, 10), collapse = ", "),
       if (n_missing > 10) paste0(" ... (", n_missing, " total)") else "")
}
phenotype_data <- phenotype_data[pheno_idx, ]

cat("Aligned", nrow(phenotype_data), "phenotype rows to PLINK sample order ",
    "(", nrow(genotypes_batch), "samples in genotype data).\n")

# ------------ Run ANOVA ----------------
ANOVA_out <- data.frame()
for (i in 1:ncol(genotypes_batch)){
  snp_data <- genotypes_batch[, i]
  # 1. Anova across groups
  anova_model <- aov(as.numeric(snp_data) ~  Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Cluster, data = phenotype_data )
  anova_summary <- summary(anova_model)
  tab <- anova_summary[[1]]
  anova_f    <- tab["Cluster", "F value"]
  anova_pval <- tab["Cluster", "Pr(>F)"]
  
  # 2. Trend test (ordinal encoding of cluster)
  trend_model <- lm(snp_data ~ numeric_cluster + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = phenotype_data)
  trend_summary <- summary(trend_model)
  trend_coef <- trend_summary$coefficients["numeric_cluster", "Estimate"]
  trend_pval <- trend_summary$coefficients["numeric_cluster", "Pr(>|t|)"]
  
  # 3. Linear regression with raw PC1_clinical
  pc1_model <- lm(snp_data ~ pc1_clinical + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = phenotype_data)
  pc1_summary <- summary(pc1_model)
  pc1_coef <- pc1_summary$coefficients["pc1_clinical", "Estimate"]
  pc1_pval <- pc1_summary$coefficients["pc1_clinical", "Pr(>|t|)"]
  
  current_res = data.frame(
    snp = snps_batch_ids[i],
    chromosome = snp_info_batch[i, "chromosome"],
    genetic_distance = snp_info_batch[i, "cM"],
    position = snp_info_batch[i, "position"],
    allele.1 = snp_info_batch[i, "allele.1"],
    allele.2 = snp_info_batch[i, "allele.2"],
    ANOVA_F = anova_f,
    ANOVA_pvalue = anova_pval,
    Trend_Estimate = trend_coef,
    Trend_pvalue = trend_pval,
    PC1_Clinical_Estimate = pc1_coef,
    PC1_Clinical_pvalue = pc1_pval
  )
  ANOVA_out = rbind(ANOVA_out, current_res)
  if (i %% 10000 == 0) {
    cat("Processed SNP", start_idx + i - 1, "\n")
    gc()
  }
}
out_file <- paste0("gwas_results_batch", batch_num, ".csv")
write.csv(ANOVA_out, out_file, row.names = FALSE)
cat("Finished writing", out_file, "\n")

rm(genetic_data_batch, genotypes_batch, ANOVA_out)
gc()
