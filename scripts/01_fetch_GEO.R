# 01_fetch_GEO.R
# Author: ychnjcpu
# Purpose: Download and prepare public IBD dataset from GEO

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"))

library(GEOquery)
library(limma)
library(tidyverse)

# ------------------------------
# Step 1: Download GEO dataset
# Example dataset: GSE48959 (IBD colonic tissue miRNA)
# You can replace this with another GEO ID
geo_id <- "GSE48959"
gset <- getGEO(geo_id, GSEMatrix = TRUE, getGPL = FALSE)

# Explicitly select the miRNA platform
gset <- gset[["GSE48959-GPL14613_series_matrix.txt.gz"]]

# ------------------------------
# Step 2: Extract expression data
exprs_data <- exprs(gset)       # expression matrix
sample_info <- pData(gset)      # sample metadata

# ------------------------------
# Step 3: Quick sanity checks
cat("Expression matrix dimensions:", dim(exprs_data), "\n")
cat("Sample metadata dimensions:", dim(sample_info), "\n")
cat("First few sample metadata rows:\n")
print(head(sample_info))

# ------------------------------
# Step 4: Save for downstream analysis
dir.create("../data", showWarnings = FALSE)  # ensure data folder exists
write.csv(exprs_data, file = "../data/expression_matrix.csv")
write.csv(sample_info, file = "../data/sample_metadata.csv")

cat("Data saved to '../data/' folder\n")