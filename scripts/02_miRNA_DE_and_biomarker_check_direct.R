# 02_miRNA_DE_and_biomarker_check.R
# Robust version — ensures proper alignment for lmFit

library(tidyverse)
library(limma)

# ------------------------------
# Step 1: Load data
exprs_data <- read.csv("../data/expression_matrix.csv", 
                       row.names = 1, 
                       check.names = FALSE, 
                       stringsAsFactors = FALSE)
sample_info <- read.csv("../data/sample_metadata.csv", stringsAsFactors = FALSE)

# Ensure column names are GSM IDs
colnames(exprs_data) <- trimws(colnames(exprs_data))
sample_info$geo_accession <- trimws(sample_info$geo_accession)
rownames(sample_info) <- sample_info$geo_accession

# ------------------------------
# Step 2: Subset only overlapping samples
common_samples <- intersect(colnames(exprs_data), rownames(sample_info))

if(length(common_samples) == 0){
  stop("No matching samples found between expression matrix and metadata!")
}

exprs_data <- exprs_data[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

# ------------------------------
# Step 3: Confirm alignment
stopifnot(
  ncol(exprs_data) == nrow(sample_info),
  all(colnames(exprs_data) == rownames(sample_info))
)

# ------------------------------
# Step 4: Define groups
group_simple <- ifelse(sample_info$`disease:ch1` == "control", "Control", "IBD")
group_simple <- factor(group_simple, levels = c("Control", "IBD"))
table(group_simple)

# ------------------------------
# Step 5: Limma differential expression
design <- model.matrix(~ 0 + group_simple)
colnames(design) <- levels(group_simple)

fit <- lmFit(exprs_data, design)
contrast.matrix <- makeContrasts(IBD_vs_Control = IBD - Control, levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

# ------------------------------
# Step 6: Save results
de_results <- topTable(fit, coef = "IBD_vs_Control", number = Inf, adjust.method = "BH")
dir.create("../results", showWarnings = FALSE)
write.csv(de_results, "../results/DE_all_miRNAs.csv")

# Extract let-7b-5p probes
let7b_rows <- grep("let-7b", rownames(de_results), ignore.case = TRUE)
let7b_results <- de_results[let7b_rows, ]
write.csv(let7b_results, "../results/DE_let7b.csv")
print(let7b_results)

# ------------------------------
# Step 7: Plot
let7b_exprs <- exprs_data[let7b_rows, , drop = FALSE]
plot_df <- data.frame(
  Expression = as.numeric(let7b_exprs),
  Group = rep(group_simple, each = nrow(let7b_exprs))
)

dir.create("../figures", showWarnings = FALSE)
png("../figures/let7b_expression_IBD_vs_Control.png", width = 700, height = 600)
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

cat("Analysis complete.\n")
cat("Results saved to:\n")
cat(" - results/DE_all_miRNAs.csv\n")
cat(" - results/DE_let7b.csv\n")
cat(" - figures/let7b_expression_IBD_vs_Control.png\n")