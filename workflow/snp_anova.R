args <- commandArgs(trailingOnly=TRUE)
plink_prefix <- args[1]
pheno_file <- args[2]
start_idx <- as.numeric(args[3])
end_idx <- as.numeric(args[4])
batch_num <- as.numeric(args[5])
# .libPaths("/path/to/your/R/library")  # Uncomment and set if using a custom R library path on your HPC
library(data.table)
library(snpStats)

phenotype_data <- fread(pheno_file)
phenotype_data$Cluster = as.factor(phenotype_data$Cluster)
phenotype_data$numeric_cluster <- as.numeric(phenotype_data$Cluster)

bim_file <- paste0(plink_prefix, ".bim")
snp_info <- fread(bim_file, header = FALSE)
colnames(snp_info) <- c("chromosome", "snp_id", "cM", "position", "allele.1", "allele.2")
rownames(snp_info) <- snp_info$snp_id

snp_info <- snp_info[!duplicated(snp_id)]
snp_ids <- snp_info$snp_id
num_snps <- length(snp_ids)

snps_batch_ids <- snp_ids[start_idx:end_idx]
snp_info_batch <- snp_info[match(snps_batch_ids, snp_info$snp_id), ] 

cat("\nReading batch", batch_num, ":", start_idx, "-", end_idx, "\n")
genetic_data_batch <- read.plink(plink_prefix, select.snps = snps_batch_ids)
genotypes_batch <- as(genetic_data_batch$genotypes, "numeric")

ANOVA_out <- data.frame()
for (i in 1:ncol(genotypes_batch)){
  snp_data <- genotypes_batch[, i]
  # 1. Anova across groups
  anova_model <- aov(as.numeric(snp_data) ~  Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Cluster, data = phenotype_data )
  anova_summary <- summary(anova_model)
  anova_f <- anova_summary[[1]][["F value"]][1]
  anova_pval <- anova_summary[[1]][["Pr(>F)"]][1]
  
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
