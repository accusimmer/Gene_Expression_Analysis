# Overview 
  This project focuses on the comparison of expression between high-impact mutations and no mutations in the 
  gene LRP1B in the primary tumor for patients that have lung squamous cell carcinoma, which could provide evidence 
  to aid in certain cancer diagnoses and treatments
  ---
# Methods
  ## Data Setup
   - RNA-seq gene expression data were retrieved from the Genomic Data Commons (GDC) using the TCGAbiolinks package in R.
   - The analysis focused on the TCGA-LUSC project with STAR - Counts workflow and included only primary tumor samples.
   - Sample metadata were filtered based on mutational status of the LRP1B gene:
     - High-impact mutation group: 31 samples
     - Non-mutated group: Randomly sampled 50 out of 115 available non-mutated samples
   - Sample IDs were matched to GDC submitter IDs to subset the expression data for both groups.
 ## Data Filtering 
   - The full expression dataset was downloaded and processed using GDCprepare(), producing a RangedSummarizedExperiment object.
   - Loci with >90% zero read counts in either group were filtered out.
   - Normalized counts were calculated using estimateSizeFactors() and extracted for variance filtering.
   - Genes were further filtered to retain only those in the top 50% of variance across samples to reduce noise in PCA.
 ## Principal Componet Analysis
   - PCA was performed using prcomp() on the variance-filtered, normalized expression matrix.
   - A visual outlier (sample ID TCGA-90-A4ED-01A-31R-A24Z-07) was removed after PCA to improve clustering.
   - The top two principal components (PC1 and PC2) were visualized with variance explained and group separation colored by mutational status using ggplot2.
##  Differential Expression Analysis
   - The filtered dataset was analyzed with DESeq2 to compare the high-impact vs. non-mutated groups.
   - Significance thresholds were defined as:
     - FDR-adjusted p-value < 0.01
     - |log₂(Fold Change)| ≥ 1.5
   - From 34,634 genes, 114 were found to be significantly differentially expressed.
##  Volcano Plot Visualization
   - A volcano plot was generated using ggplot2 and ggrepel, with:
     - X-axis: log₂(fold change)
     - Y-axis: –log₁₀(FDR-adjusted p-value)
     - Color scale based on mean gene expression
     - Top 10 genes with smallest p-values labeled
     - Threshold lines:
        - Vertical: log₂FC ±1.5 
        - Horizontal: p = 0.05 








