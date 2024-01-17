# Colon-Adenocarcinoma
This repository contains matadata and codes necessary for analysis of the scRNA-seq data of human Small Intestine and Colon presented in the artical titled Heterogeneous Immune Landscape between Small Intestine and Colon Provides Novel Insights of Intraepithelial T Cells in the Progression of Colon Adenocarcinoma.

## Data processing

### Single Cell RNASeq
For scRNASeq data, a multi-step filtering process was further performed to remove likely cell debris and doublets during the process of unsupervised clustering. Clusters with nonspecific markers, significantly low number of expressed genes and high percentage of transcripts from the mitochondrial genome were considered low-quality cells or likely debris and thus excluded. We then used the DoubletFinder (v2.0.3) algorithm to calculate the doublet score for each cell. Clusters with high doublet score, expressing hybrid cell lineage markers and exhibiting an aberrantly high gene count were considered likely doublets and thus excluded. To further clean doublets that could have been missed in the above steps, we generated UMAP plots and carefully reviewed expression of canonical marker genes across cell clusters, and cells co-expressing markers of discrepant lineages or cells matched with TCR/BCR annotations but expressing discrepant markers were further excluded. After these strict multi-step filtering process, of the 534,612 cells that passed the initial quality filtering, 22,104 cells were further excluded for low-quality, 45,786 cells were further excluded for possible doublets, and 802 enterocytes were also excluded, resulting in 465,920 immune cells for down-stream analyses.

### Bulk RNASeq
To prepare data for downstream analysis, FastQC (v0.12.1) was used to assess the quality of raw sequencing data, and fastp software (v0.23.4) was applied to remove low-quality reads and any residual adapter sequences. The high-quality reads were aligned to human genome reference GRCh38 (hg38) using the STAR algorithm (v2.7.10b). ResQC (v5.0.3) is applied to assess the quality of alignment results using various alignment metrics, including coverage, mapping quality, and distribution of mapped reads. Then, FeatureCount (v2.0.6) was utilized to extract the feature count matrix from alignments. 

## Downloading the data
Cell level metadata is available in the provided /sc_data/meta.Rdata, which contains sample, tissue type, major cell types and detailed cell types.

Bulk level metadata is available in the provided /bulk_data/meta.Rdata, which contains sample and tissue type.

The raw data arw uploaded to GSA: scRNASeq: HRA006401, bulk RNASeq: HRA006350

The processed data are uploaded to mendeley: scRNASeq: 10.17632/6czch25jyb.1, bulk RNASeq: 10.17632/hb9jjk2gbz.1.
