This README explains the steps and code used to analyze RNA-Seq data, perform differential gene expression analysis, and generate survival models for breast cancer patients with ERBB2 amplification.

Overview
The pipeline processes RNA-Seq, copy number aberration (CNA), and clinical patient data. It identifies differentially expressed genes between ERBB2-amplified and non-amplified patients, performs pathway enrichment analysis, visualizes results, and builds survival models.

Prerequisites
Ensure the following R packages are installed:

Data manipulation: BiocManager, DESeq2, ggplot2, pheatmap
Pathway analysis: clusterProfiler, org.Hs.eg.db, pathview
Survival analysis: glmnet, survival

Steps in the Analysis

1. Dataset Preparation
The dataset (brca_tcga_pan_can_atlas_2018.tar.gz) is extracted to the working directory.
The following files are loaded:
data_mrna_seq_v2_rsem.txt (RNA-Seq data)
data_clinical_patient.txt (patient metadata)
data_cna.txt (CNA data)

2. Matching Patient IDs
RNA-Seq, CNA, and clinical data are filtered to include only common patient IDs for downstream analysis.

3. Metadata Creation
Metadata is generated based on ERBB2 CNA levels. Patients with CNA levels > 0 for the ERBB2 gene are labeled as amplified.

4. Normalization and Differential Expression Analysis
Normalization and differential expression analysis are performed using DESeq2.
Results include:
Top 10 differentially expressed genes ranked by adjusted p-value.
Top 10 genes ranked by log2 fold change.

5. Pathway Enrichment Analysis
Significant genes (adjusted p-value < 0.05) are mapped to KEGG pathways using clusterProfiler.
Pathway enrichment results include:
Bar plots of top enriched pathways.
KEGG pathway IDs and their significance levels.

6. Variance-Stabilized Transformation and Visualization
Variance-stabilized transformed (VST) values are computed for PCA and heatmap generation:
PCA plot: Visualizes clustering of samples by ERBB2 amplification.
Heatmap: Displays the top 10 variable genes across samples.

7. Survival Analysis
Survival models are built using the glmnet package:
Survival data (OS_MONTHS and OS_STATUS) is used to compute survival probabilities.
High- and low-risk groups are defined based on median risk scores.
Kaplan-Meier survival curves are plotted to compare survival probabilities.

Outputs

Differential Expression Results

Top significant genes ranked by adjusted p-value and fold change.
Pathway Enrichment

Bar plots of top enriched KEGG pathways.

PCA plot grouped by ERBB2 amplification.

Heatmap of the top 10 variable genes.

Kaplan-Meier survival curve for high- and low-risk groups.

How to Run the Code

Set the correct file paths for the dataset in the downloads_path variable.

Run the code in R studio or another R environment with all required packages installed.

Ensure internet access for package installation and KEGG pathway queries.

Key Notes

ERBB2 Amplification: Used as the grouping variable for differential expression analysis and survival models.

Log2 Fold Change: Reflects the magnitude of expression differences between amplified and non-amplified groups.

Pathway Enrichment: Highlights biological processes associated with differentially expressed genes.

Kaplan-Meier Analysis: Assesses survival differences based on high- and low-risk scores from gene expression data.
