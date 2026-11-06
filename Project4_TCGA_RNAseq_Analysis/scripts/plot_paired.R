
#!/usr/bin/env Rscript

# Paired NKX2-1 expression visualization
# TCGA-LUAD: Primary Tumor vs Solid Tissue Normal

library(ggplot2)

input_file <- "results/NKX2-1_paired_expression.tsv"

d <- read.delim(input_file, stringsAsFactors = FALSE)

# Convert expression values
d$tumour_log2 <- as.numeric(d$tumour_log2)
d$normal_log2 <- as.numeric(d$normal_log2)

# Check data
if (nrow(d) == 0) {
  stop("No paired samples available.")
}

if (anyNA(d$tumour_log2) || anyNA(d$normal_log2)) {
  stop("Missing or invalid expression values.")
}

# Convert paired data into long format
plot_data <- rbind(
  data.frame(
    case_id = d$case_id,
    sample_type = "Solid Tissue Normal",
    expression = d$normal_log2
  ),
  data.frame(
    case_id = d$case_id,
    sample_type = "Primary Tumor",
    expression = d$tumour_log2
  )
)

plot_data$sample_type <- factor(
  plot_data$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)

# Paired Wilcoxon test
wt <- wilcox.test(
  d$tumour_log2,
  d$normal_log2,
  paired = TRUE,
  exact = FALSE
)

n <- nrow(d)
y_max <- max(plot_data$expression)

# Create paired boxplot
p <- ggplot(
  plot_data,
  aes(x = sample_type, y = expression)
) +

  geom_line(
    aes(group = case_id),
    color = "grey65",
    alpha = 0.6
  ) +

  geom_boxplot(
    width = 0.35,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +

  geom_point(
    position = position_jitter(width = 0.06, height = 0),
    size = 1.8,
    alpha = 0.65
  ) +

  annotate(
    "text",
    x = 1:2,
    y = y_max + 0.4,
    label = paste0("n = ", n)
  ) +

  labs(
    title = "NKX2-1 Expression in TCGA-LUAD",
    subtitle = sprintf(
      "Paired Wilcoxon signed-rank test: p = %.4f",
      wt$p.value
    ),
    x = NULL,
    y = "Expression, log2(TPM + 1)"
  ) +

  coord_cartesian(
    ylim = c(min(plot_data$expression), y_max + 1)
  ) +

  theme_classic(base_size = 13) +

  theme(
    axis.text.x = element_text(color = "black"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# Save figures
ggsave(
  "results/NKX2-1_paired_boxplot.pdf",
  plot = p,
  width = 6,
  height = 5
)

ggsave(
  "results/NKX2-1_paired_boxplot.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

print(p)

cat("Paired plot saved successfully.\n")
