
#!/usr/bin/env Rscript

# Paired NKX2-1 expression analysis
# TCGA-LUAD: Primary Tumor vs Solid Tissue Normal

input_file <- "results/NKX2-1_paired_expression.tsv"

# Read paired expression data
d <- read.delim(input_file, stringsAsFactors = FALSE)

# Check required columns
required <- c(
  "case_id",
  "tumour_tpm",
  "normal_tpm",
  "tumour_log2",
  "normal_log2"
)

if (!all(required %in% names(d))) {
  stop("Input file is missing required columns.")
}

if (nrow(d) == 0) {
  stop("No paired samples found.")
}

# Convert expression columns to numeric
d$tumour_tpm <- as.numeric(d$tumour_tpm)
d$normal_tpm <- as.numeric(d$normal_tpm)
d$tumour_log2 <- as.numeric(d$tumour_log2)
d$normal_log2 <- as.numeric(d$normal_log2)

# Remove incomplete pairs
d <- d[complete.cases(d[, required]), ]

# Paired Wilcoxon signed-rank test
wt <- wilcox.test(
  d$tumour_log2,
  d$normal_log2,
  paired = TRUE,
  exact = FALSE,
  conf.int = TRUE
)

# Calculate paired differences
differences <- d$tumour_log2 - d$normal_log2

# Calculate medians
tumour_median <- median(d$tumour_tpm)
normal_median <- median(d$normal_tpm)

cat("NKX2-1: Paired Expression Analysis\n")
cat("----------------------------------\n")

cat("Matched patients:", nrow(d), "\n")

cat("\nMedian TPM:\n")
cat("Tumour:", tumour_median, "\n")
cat("Normal:", normal_median, "\n")

cat("\nMedian paired log2 difference:",
    median(differences), "\n")

cat("Paired fold-change estimate (log2 scale):",
    2^median(differences), "\n")

cat("\nWilcoxon signed-rank test:\n")
cat("P-value:",
    format.pval(wt$p.value, digits = 5), "\n")

cat("\nHodges-Lehmann location-shift estimate:\n")
print(wt$estimate)

cat("\n95% confidence interval:\n")
print(wt$conf.int)
