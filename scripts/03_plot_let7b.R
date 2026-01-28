# 03_plot_let7b.R
# Boxplot of let-7b expression in IBD vs Control
# Requires: exprs_data and group_simple already in memory

library(ggplot2)

# ------------------------------
# Ensure the figures folder exists
dir.create("C:/projects/IBD_miRNA_biomarker/figures", showWarnings = FALSE)

# ------------------------------
# Extract let-7b expression
let7b_rows <- grep("let-7b", rownames(exprs_data), ignore.case = TRUE)
let7b_exprs <- exprs_data[let7b_rows, , drop = FALSE]

# Prepare data frame for plotting
plot_df <- data.frame(
  Expression = as.numeric(let7b_exprs),
  Group = rep(group_simple, each = nrow(let7b_exprs))
)

# ------------------------------
# Plot and save
png("C:/projects/IBD_miRNA_biomarker/figures/let7b_expression_IBD_vs_Control.png",
    width = 700, height = 600)
ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "let-7b expression in IBD vs Control",
    subtitle = "GSE48959 (miRNA, colon tissue)",
    y = "Normalized expression"
  ) +
  scale_fill_brewer(palette = "Set2")
dev.off()

cat("Boxplot saved to: C:/projects/IBD_miRNA_biomarker/figures/let7b_expression_IBD_vs_Control.png\n")
