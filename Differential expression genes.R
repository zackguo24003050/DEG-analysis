# Load libraries
library(readr)
library(limma)

# Read clinical data
pheno <- read_tsv("/Users/guozekai/Downloads/pheno_sclc_ucologne_2015.tsv")

# Stage classification
pheno <- pheno[!is.na(pheno$uicc_tumor_stage), ]
pheno$stage_group <- ifelse(
  pheno$uicc_tumor_stage %in% c("I", "Ia", "Ib", "II", "IIa", "IIb"), "early",
  ifelse(pheno$uicc_tumor_stage %in% c("III", "IIIa", "IIIb", "IV"), "advanced", NA)
)
pheno <- pheno[!is.na(pheno$stage_group), ]

# Read expression data
expr <- readRDS("/Users/guozekai/expr_normalized_median_polish.rds")

# Match clinical and expression samples
common_samples <- intersect(colnames(expr), pheno$patient_id)
expr <- expr[, common_samples]
pheno <- pheno[match(common_samples, pheno$patient_id), ]

# Create design matrix
group <- factor(pheno$stage_group, levels = c("early", "advanced"))
design <- model.matrix(~ group)

# Differential expression analysis
fit <- lmFit(expr, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "groupadvanced", number = Inf, sort.by = "P")

# Filter DEG
deg_explore <- results[results$P.Value < 0.01 & abs(results$logFC) > 1, ]
cat("Exploratory DEGs (p < 0.01 & |logFC| > 1):", nrow(deg_explore), "\n")

# Package for plotting
library(ggplot2)

# Blue color for down-regulated and red color for up-regulated
results$threshold <- "Not Sig"
results$threshold[results$P.Value < 0.05 & results$logFC > 1] <- "Up"
results$threshold[results$P.Value < 0.05 & results$logFC < -1] <- "Down"

# Y-axis is -log10(p-value)
results$logP <- -log10(results$P.Value)

# Plot graph
ggplot(results, aes(x = logFC, y = logP, color = threshold)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value"
  ) +
  theme_minimal()


