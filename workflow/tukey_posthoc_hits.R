args <- commandArgs(trailingOnly = TRUE)

plink_prefix <- args[1]
pheno_file <- args[2]
hits_file <- args[3]
out_file <- args[4]

library(data.table)
library(snpStats)

phenotype_data <- fread(pheno_file)
phenotype_data$Cluster <- as.factor(phenotype_data$Cluster)

hits <- fread(hits_file)

if (!"snp" %in% colnames(hits)) {
  stop("hits_file must contain a column named 'snp'")
}

hit_snps <- unique(hits$snp)

cat("Reading", length(hit_snps), "significant SNPs for Tukey post-hoc testing\n")

genetic_data <- read.plink(plink_prefix, select.snps = hit_snps)
genotypes <- as(genetic_data$genotypes, "numeric")

tukey_out <- data.frame()

for (i in seq_along(hit_snps)) {
  snp_id <- hit_snps[i]
  snp_data <- genotypes[, i]

  anova_model <- aov(
    as.numeric(snp_data) ~ Age + Sex +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
      Cluster,
    data = phenotype_data
  )

  tukey_res <- tryCatch(
    TukeyHSD(anova_model, which = "Cluster"),
    error = function(e) NULL
  )

  if (!is.null(tukey_res) && "Cluster" %in% names(tukey_res)) {
    tukey_tab <- as.data.frame(tukey_res$Cluster)
    tukey_tab$contrast <- rownames(tukey_tab)

    current <- data.frame(
      snp = snp_id,
      contrast = tukey_tab$contrast,
      diff = tukey_tab$diff,
      lwr = tukey_tab$lwr,
      upr = tukey_tab$upr,
      tukey_adj_pvalue = tukey_tab[["p adj"]]
    )

    tukey_out <- rbind(tukey_out, current)
  }

  if (i %% 100 == 0) {
    cat("Processed Tukey post-hoc SNP", i, "of", length(hit_snps), "\n")
    gc()
  }
}

write.csv(tukey_out, out_file, row.names = FALSE)

cat("Finished writing", out_file, "\n")