# Gene Differential Expression Analysis

An R coursework project exploring gene expression in the TCGA breast cancer dataset. The analysis compares patients grouped by ERBB2 copy number, looks at pathway enrichment, and explores survival patterns.

## What's here?

| File | Purpose |
| --- | --- |
| [Analysis script](scripts/Bio_Assignment_Fin.R) | Data preparation, differential expression, pathway enrichment, plots, and survival modelling. |
| [Maths pre-reading notebook](references/Mathematics_for_Neuroscience_pre_read.ipynb) | Supplementary neuroscience maths material; separate from the R analysis. |

The dataset and generated results are not included.

## Get started

You will need R and an R environment such as RStudio. The script uses these packages:

- Bioconductor: `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, and `pathview`.
- CRAN: `BiocManager`, `ggplot2`, `pheatmap`, `glmnet`, and `survival`.

Install them from the R console:

```r
install.packages(c("BiocManager", "ggplot2", "pheatmap", "glmnet", "survival"))
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "pathview"))
```

Some sections of the script also call package installers. Internet access is needed for installation and KEGG queries.

### Prepare the data

Obtain the `brca_tcga_pan_can_atlas_2018.tar.gz` archive separately. Put it in the repository's root folder, or set `downloads_path` near the top of the script to the folder containing your archive.

The script expects the archive to extract into a folder named `brca_tcga_pan_can_atlas_2018` containing:

```text
data_mrna_seq_v2_rsem.txt
data_clinical_patient.txt
data_cna.txt
```

### Run the analysis

1. Open [scripts/Bio_Assignment_Fin.R](scripts/Bio_Assignment_Fin.R) in RStudio.
2. Set your working directory to the repository's root folder.
3. Check `downloads_path`, then run the numbered sections in order.
4. Inspect the intermediate tables and plots as you go.

The script switches the working directory to the extracted data folder. Reset it to the repository root before running the script again.

## Analysis workflow

1. Load the expression, copy number, and clinical tables.
2. Match their patient IDs and combine duplicate gene-symbol rows.
3. Label patients using ERBB2 copy number. The script calls any value greater than zero "amplified".
4. Prepare and round the expression values, then run DESeq2.
5. Rank genes by adjusted p-value and absolute log2 fold change.
6. Use genes with adjusted p-value below 0.05 for KEGG enrichment.
7. Create a PCA plot and a heatmap from variance-stabilized expression values.
8. Fit a cross-validated penalized Cox model and plot survival for groups split at the median fitted risk score.

## What to expect

Results are shown in the R console and plot viewer: gene rankings, pathway summaries, a pathway bar plot, PCA, a heatmap, and Kaplan-Meier curves. The script does not currently export these automatically; save the tables and plots you want to keep.

## Current limitations

This is an exploratory coursework script. It has no bundled dataset, automated tests, or recorded package versions.

Before interpreting results, check that sample ordering agrees across the expression, copy number, and clinical tables. The script rounds the provided expression values and uses a positive-CNA grouping rule; these choices need review for any further analysis. The survival curves use risk scores fitted on the same patients and do not measure performance on a separate test cohort.

## Supplementary notebook

The pre-reading notebook covers differential equations and linear algebra. Its header credits Barry Dillon and notes by Áine Byrne and Rachel Nicks. It uses Python, NumPy, and Matplotlib.

The image-matrix exercise loads `pixel_matrix.npy`, which is not included. Obtain that file separately and place it alongside the notebook, or skip that exercise. Some cells deliberately demonstrate invalid matrix operations and will raise an error.
