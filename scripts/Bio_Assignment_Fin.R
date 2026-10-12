# Breast cancer gene expression analysis
# Run the sections in order in RStudio so you can inspect the tables and plots.
# Put the dataset archive in your working directory, or change downloads_path below.

# 1. Choose the folder containing the downloaded dataset archive.
downloads_path = getwd()
file_path = paste(downloads_path,"brca_tcga_pan_can_atlas_2018.tar.gz", sep = "/" )

# 2. Extract the archive and work inside its data folder.
untar(file_path)
folder_path = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/")
setwd(new_dir)

# 3. Load the gene expression table.

data_Rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# 4. Load patient details, skipping the file's metadata comments.

data_patient <- read.delim("data_clinical_patient.txt",
                           stringsAsFactors = FALSE,
                           header = TRUE,
                           comment.char = "#")

# 5. Load the copy number table.

data_cna = read.delim("data_cna.txt", stringsAsFactors = FALSE, header = TRUE)

# 6. Keep patients present in all three tables.

# Remove rows with no Hugo_Symbol
rna_filtered <- data_Rnaseq[data_Rnaseq$Hugo_Symbol != "", ]
cna_filtered <- data_cna[data_cna$Hugo_Symbol != "", ]

# Combine rows that share a gene symbol by summing their numeric values.
rna_filtered <- aggregate(. ~ Hugo_Symbol, data = rna_filtered, FUN = sum)
cna_filtered <- aggregate(. ~ Hugo_Symbol, data = cna_filtered, FUN = sum)

# R changes hyphens to dots in column names; restore the patient barcode format.
rna_sample_ids <- substr(gsub("\\.", "-", colnames(rna_filtered)[3:ncol(rna_filtered)]), 1, 12)
cna_sample_ids <- substr(gsub("\\.", "-", colnames(cna_filtered)[3:ncol(cna_filtered)]), 1, 12)

# Use the clinical table's patient IDs as the matching reference.
patient_ids <- data_patient$PATIENT_ID

# Only keep patients with expression, copy number, and clinical data.
common_ids <- intersect(intersect(rna_sample_ids, cna_sample_ids), patient_ids)

# Apply the shared patient list to each table.
rna_filtered <- rna_filtered[, c("Hugo_Symbol", "Entrez_Gene_Id", colnames(rna_filtered)[substr(gsub("\\.", "-", colnames(rna_filtered)), 1, 12) %in% common_ids])]
cna_filtered <- cna_filtered[, c("Hugo_Symbol", "Entrez_Gene_Id", colnames(cna_filtered)[substr(gsub("\\.", "-", colnames(cna_filtered)), 1, 12) %in% common_ids])]
patient_filtered <- data_patient[data_patient$PATIENT_ID %in% common_ids, ]

# 7. Build the ERBB2 grouping labels from copy number values.

metadata <- matrix(0, dim(cna_filtered)[2], 1)

pat_ids <- data_patient$PATIENT_ID

for (i in 1:dim(cna_filtered)[2]) {
  pat_barcode <- colnames(cna_filtered)[i]
  pat_barcode <- substr(gsub("\\.", "-", pat_barcode), 1, 12)
  idx <- which(pat_barcode == pat_ids)
  if (length(idx) > 0) {
    erbb2_cna_value <- as.numeric(cna_filtered[cna_filtered$Hugo_Symbol == "ERBB2", i])
    metadata[i, 1] <- as.numeric(erbb2_cna_value > 0)  # This script uses any positive CNA value for the amplified group.
  }
}

colnames(metadata) <- "ERBB2_Amplified"

metadata[is.na(metadata)] =0

# 8. Prepare the expression matrix for DESeq2.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DESeq2 before running the expression analysis.

BiocManager::install("DESeq2")

library(DESeq2)

# Drop the gene-name and Entrez-ID columns from the expression matrix.
assay <- rna_filtered[, -(1:2)]
rownames(assay) <- rna_filtered$Hugo_Symbol

# Ensure the data is numeric
assay <- apply(assay, 2, as.numeric)
rownames(assay) <- rna_filtered$Hugo_Symbol
assay <- as.matrix(assay)

# Replace NA values with 0 and set negative values to 0
assay[is.na(assay)] <- 0
assay[assay < 0] <- 0

# Filter out genes with too few counts
smallestGroupSize <- 3
keep <- rowSums(assay >= 10) >= smallestGroupSize
assay <- assay[keep, ]

# Remove the two labels created for the non-sample columns.
assay_sample_ids <- colnames(assay)

metadata <- metadata[-c(1, 2), , drop = FALSE]

# 9. Compare gene expression between the two ERBB2 groups.

assay <- round(assay)

metadata <- as.data.frame(metadata)

# Convert ERBB2_Amplified to a factor
metadata$ERBB2_Amplified <- factor(metadata$ERBB2_Amplified, levels = c(0, 1))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = metadata,
                              design = ~ ERBB2_Amplified)

# Normalize data
dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

res = results(dds)

# Show the ten genes with the smallest adjusted p-values.

res[order(res$padj)[1:10],]

top_genes_by_fold_change <- res[order(abs(res$log2FoldChange), decreasing = TRUE)[1:10], ]

# Print the top 10 genes
top_genes_by_fold_change

# Get the total number of genes analyzed
total_genes <- nrow(dds)
cat("Total number of genes analyzed:", total_genes, "\n")

# 10. Look for KEGG pathways represented among significant genes.

# Install required packages if not already installed

install.packages("BiocManager", dependencies = TRUE)
library(BiocManager)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Extract significant genes based on adjusted p-value cutoff (e.g., < 0.05)
significant_genes <- rownames(res[which(res$padj < 0.05), ])

# Map Hugo_Symbols to Entrez IDs
gene_entrez_ids <- bitr(significant_genes,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

# Keep each Entrez ID once so a gene is not counted twice.
gene_entrez_ids <- gene_entrez_ids[!duplicated(gene_entrez_ids$ENTREZID), ]

# Test for enrichment against human KEGG pathways.
kegg_results <- enrichKEGG(gene = gene_entrez_ids$ENTREZID,
                           organism = "hsa",
                           pvalueCutoff = 0.05)

# Display top 10 enriched pathways
head(kegg_results, n = 10)

# Bar plot of enriched KEGG pathways
barplot(kegg_results, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")

summary(res$padj)

# 11. Transform expression values for exploratory plots.

# Variance-stabilizing transformation (VST)
vsd <- vst(dds, blind = TRUE)  # 'blind = TRUE' for unsupervised transformation

# Get the transformed expression values
vst_counts <- assay(vsd)  # Extract VST-transformed values as a matrix

# 12. Plot sample similarity and the most variable genes.

library(ggplot2)
library(pheatmap)

par(mfrow = c(1, 1))

# PCA plot using ERBB2 amplified
pca_data <- plotPCA(vsd, intgroup = "ERBB2_Amplified", returnData = TRUE)

ggplot(pca_data, aes(x = PC1, y = PC2, color = ERBB2_Amplified)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Plot Grouped by ERBB2 Amplification",
    x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
    y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# Top variable genes
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 10)

# Get top variable genes from VST
heatmap_data <- assay(vsd)[top_var_genes, ]

# Label each heatmap column with its ERBB2 group.
annotation_col <- as.data.frame(colData(vsd)[, "ERBB2_Amplified", drop = FALSE])
rownames(annotation_col) <- colnames(vsd)

# Create heatmap
pheatmap(
  heatmap_data,
  scale = "row",  # Scale each gene (row) to Z-scores
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "Heatmap of Top 10 Variable Genes"
)

# 13. Fit a survival model and compare the resulting risk groups.

install.packages("glmnet")
install.packages("survival")

library(glmnet)
library(survival)

# Read follow-up time and turn the recorded status into an event flag.
survival_time <- as.numeric(patient_filtered$OS_MONTHS)
survival_status <- as.numeric(patient_filtered$OS_STATUS == "1:DECEASED")

# Filter VST values for differentially expressed genes
de_genes <- rownames(res[res$padj < 0.05 & !is.na(res$padj), ])
vst_de_genes <- assay(vsd)[de_genes, ]
colnames(vst_de_genes) <- substr(gsub("\\.", "-", colnames(vst_de_genes)), 1, 12)

# Find common samples
common_samples <- intersect(colnames(vst_de_genes), patient_filtered$PATIENT_ID)
vst_de_genes <- vst_de_genes[, common_samples]
survival_time <- survival_time[match(common_samples, patient_filtered$PATIENT_ID)]
survival_status <- survival_status[match(common_samples, patient_filtered$PATIENT_ID)]

# glmnet expects patients in rows and genes in columns.
predictor_matrix <- t(vst_de_genes)

# Swap zero or negative times with a small positive value
survival_time[survival_time <= 0] <- 0.01

valid_samples <- survival_time > 0
survival_time <- survival_time[valid_samples]
survival_status <- survival_status[valid_samples]
predictor_matrix <- predictor_matrix[valid_samples, ]

# Combine follow-up time and event status for the Cox model.
response <- Surv(survival_time, survival_status)

# Fit glmnet with cross-validation to find optimal penalty parameter (lambda)
cv_fit <- cv.glmnet(
  x = predictor_matrix,
  y = response,
  family = "cox",
  alpha = 1
)

# Optimal lambda
optimal_lambda <- cv_fit$lambda.min
cat("Optimal lambda:", optimal_lambda, "\n")

# Refit the model using the penalty selected by cross-validation.
final_model <- glmnet(
  x = predictor_matrix,
  y = response,
  family = "cox",
  alpha = 1,
  lambda = optimal_lambda
)

# Predict risk scores
risk_scores <- predict(final_model, newx = predictor_matrix, type = "link")

# Kaplan-Meier survival plot
high_risk <- risk_scores > median(risk_scores)
surv_fit <- survfit(Surv(survival_time, survival_status) ~ high_risk)

plot(surv_fit, col = c("blue", "red"), lwd = 2,
     xlab = "Time (months)", ylab = "Survival Probability",
     main = "Kaplan-Meier Survival by Risk Group")
legend("bottomleft", legend = c("Low Risk", "High Risk"), col = c("blue", "red"), lwd = 2)
