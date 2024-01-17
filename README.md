# Colon-Adenocarcinoma
This repository contains some data and codes necessary for analysis of the scRNA-seq data of human Small Intestine and Colon presented in the artical titled _"Heterogeneous Immune Landscape between Small Intestine and Colon Provides Novel Insights of Intraepithelial T Cells in the Progression of Colon Adenocarcinoma"_.

## Data processing

### Single Cell Data
Raw scRNA-seq data were preprocessed using CellRanger Single Cell Software Suite (v6.1.2, 10X Genomics) by maping to human reference genome build 38. Cells expressing fewer than 200 genes or more than 5,000 genes were considered low-quality or potential doublets and were removed. Likely apoptotic cells where >10% transcripts were derived the mitochondrial genome were also excluded. scTCR-seq/scBCR-seq data were preprocessed by CellRanger v6.1.2 for V(D)J sequence assembly and TCR/BCR reconstruction using the human reference genome build 38. Only high-confidence and productive TRA/TRB annotations were used for further analysis. After that, DoubletFinder (v2.0.3) algorithm was further used to exclude doublets.
### Bulk RNASeq
To prepare data for downstream analysis, FastQC (v0.12.1) was used to assess the quality of raw sequencing data, and fastp software (v0.23.4) was applied to remove low-quality reads and any residual adapter sequences. The high-quality reads were aligned to human genome reference GRCh38 (hg38) using the STAR algorithm (v2.7.10b). ResQC (v5.0.3) is applied to assess the quality of alignment results using various alignment metrics, including coverage, mapping quality, and distribution of mapped reads. Then, FeatureCount (v2.0.6) was utilized to extract the feature count matrix from alignments. 

## Downloading the data
Specific input data for each figure is also inclued in /input_data.

The raw data are uploaded to GSA: scRNASeq: HRA006401, bulk RNASeq: HRA006350  
The processed data are uploaded to mendeley: scRNASeq: 10.17632/6czch25jyb.1, bulk RNASeq: 10.17632/hb9jjk2gbz.1.

## Data visualization
### Requirements
Tested on macOS Big Sur  
1. R version 4.3.1
2. R packages
   - Seurat
   - ggplot2
   - data.table
   - dplyr
   - tidyr
   - ArchR
   - pheatmap
   - ggpubr
