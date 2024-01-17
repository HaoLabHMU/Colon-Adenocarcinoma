library(dplyr);library(ggplot2);library(pheatmap)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer);library(ggthemes)
#----------------------------------------------------------------------------------------------------------
####  load functions  ####
#----------------------------------------------------------------------------------------------------------
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
               colorRampPalette(c("#fddbc7","#b2182b"))(84))
bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
groupMeans <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}
# TCGA and Gtex datasets were obtained as described in the manuscript
#----------------------------------------------------------------------------------------------------------
#### Clustering analysis (Fig.1A, 1B) ####
#----------------------------------------------------------------------------------------------------------
load("TCGA_Gtex.rda")#merged gene exp, TCGA CRC clinical information and definition of Gtex samples
library(tidyverse); library(ggdendro); library(RColorBrewer); library(plotly)
library(dplyr);library(dendextend)
Tissue_color <- c(GTEx_SI = "#1a9850", GTEx_Sigmoid = "#74add1", GTEx_Transverse="#8073ac",
                   TCGA_CA = "#d6604d",TCGA_colon="#fddbc7")

# top 3000 CD45 correlated genes
{
  CD45 <- TG_exp["PTPRC",]
  CD45.cor <- cor(CD45,t(TG_exp),use = "na.or.complete")
  CD45.cor <- CD45.cor[1,]
  CD45.cor <- CD45.cor[!is.na(CD45.cor)]
  CD45.cor <- CD45.cor[order(CD45.cor,decreasing = T)]
  imm.genes <- names(CD45.cor)[1:3000]#[CD45.cor>0.4]
}

dist_mat <- dist(t(TG_exp[imm.genes,]), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'ward.D2')
avg_dend_obj <- as.dendrogram(hclust_avg)
dendrogram_data <- dendro_data(avg_dend_obj)
species_color <- Tissue_color[TG_Samples[dendrogram_data$labels$label,"SampleType"]]
avg_dend_obj %>% set("leaves_pch", c(0)) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", species_color) %>% #node point color
  plot(main = "Leaves points")

# survival differences of CRCs
{
  avg_dend_obj %>% color_branches(.,k = 2,groupLabels = F) %>% plot
  desired_branch <- dendextend::cutree(
    avg_dend_obj,
    k = 2,
    order_clusters_as_data = FALSE)
  
  library(phenoTest);library(survival);library(survminer);
  CRC_clin$group <- desired_branch[rownames(CRC_clin)]
  fit <- survfit(Surv(DFI.time, DFI) ~ group,data = CRC_clin)
  ggsurvplot(fit, data = CRC_clin, risk.table = TRUE,pval = TRUE,risk.table.y.text = FALSE,
             xscale=365,break.x.by=365*2,xlim=c(0,365*13))
}
#----------------------------------------------------------------------------------------------------------
#### Heatmap of IEL markers (Figure 1F) ####
#----------------------------------------------------------------------------------------------------------
IEL.markers <- c("CD160","TMIGD2","IKZF2","ENTPD1","ID3",
                 "ITGAE","KLRC1","KLRC2","KLRC3","KIR2DL4")
df <- TG_Samples
df$clusters <- NA
df$clusters[df$SampleType %in% c("GTEx_Transverse","TCGA_colon","GTEx_Sigmoid")] <- "Colon"
df$clusters[df$SampleType %in% c("GTEx_SI")] <- "GTEx_SI"
df$clusters[df$group==2 & df$SampleType=="TCGA_CA"] <- "CA_C2"
df$clusters[df$group==1 & df$SampleType=="TCGA_CA"] <- "CA_C1"
df$clusters[df$group==1 & df$clusters=="Colon"] <- "Colon_C1"
df$clusters[df$group==2 & df$clusters=="Colon"] <- "Colon_C2"
df <- df[order(df$clusters),]

expMat <- TG_exp[IEL.markers,rownames(df)]
pheatmap(expMat, cluster_rows=F, cluster_cols=F,
         scale="row",show_colnames=F,show_rownames=T,
         breaks = bk,color = colour_bk,angle_col=45,
         annotation_col = df[,"clusters",drop=F])

#----------------------------------------------------------------------------------------------------------
#### Figure 1C ####
#----------------------------------------------------------------------------------------------------------
#1 DEGs from CA_C2 vs CA_C1 and ileum vs colon:
{
  CRC_clin$group <- desired_branch[rownames(CRC_clin)]
  CRC_clin <- CRC_clin[order(CRC_clin$group),]
  expMat <- TG_exp[,rownames(CRC_clin)]
  expMat <- groupMeans(expMat,groups=CRC_clin$group,sparse = F)
  df <- data.frame(rownames(expMat));colnames(df) <- "gene"
  rownames(df) <- df$gene
  df$FC_CA <- expMat[,2]-expMat[,1]
  
  df2 <- TG_Samples[TG_Samples$SampleType != "TCGA_CA",]
  df2$SampleType[df2$SampleType %in% c("GTEx_Transverse","TCGA_colon","GTEx_Sigmoid")] <- "Colon"
  expMat <- TG_exp[,rownames(df2)]
  expMat <- groupMeans(expMat,groups=df2$SampleType,sparse = F)
  df$FC_IC <- expMat[,"GTEx_SI"] - expMat[,"Colon"]
  
  # only keep protein coding genes
  iOrd <- intersect(geneTable$gene_name[geneTable$gene_type=="protein_coding"],df$gene)
  df <- df[iOrd,]
}

#2. immune genes
{
  df$label <- NA
  iOrd1 <- which(df$FC_CA>0.5 & df$FC_IC>1) 
  iOrd2 <- which(df$FC_CA< -0.5 & df$FC_IC < -1)
  df$label[iOrd1] <- df$gene[iOrd1]
  
  #plot for immune genes:
  library(ggrepel)
  df2 <- df[df$gene %in% imm.genes,]
  cor.test(df2$FC_CA,df2$FC_IC)
  p <- ggplot(data=df2, aes(x = FC_IC, y = FC_CA)) +
    geom_point(stroke = 0, alpha = 0.8,shape = 16,size=3) + 
    theme_classic() +
    #geom_text_repel(colour="black") +
    scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_vline(xintercept=0, linetype="dashed", color = "red")+
    xlim(-0.5,5)+ylim(-0.5,2)
}

#3. other genes
{
  df2 <- df[df$gene %ni% imm.genes,]
  cor.test(df2$FC_CA,df2$FC_IC)
  p <- ggplot(data=df2, aes(x = FC_IC, y = FC_CA)) +
    geom_point(stroke = 0, alpha = 0.8,shape = 16,size=3) + 
    theme_classic() +
    #geom_text_repel(colour="black") +
    scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_vline(xintercept=0, linetype="dashed", color = "red")+
    xlim(-3.2,5)+ylim(-1,2)
}
#----------------------------------------------------------------------------------------------------------
####signals using ssGSEA (Fig.1E,1G, 1I) ####
#----------------------------------------------------------------------------------------------------------
IEL.markers <- c("CD160","TMIGD2","IKZF2","ENTPD1","ID3",
                 "ITGAE","KLRC1","KLRC2","KLRC3","KIR2DL4")
T_act <- fread("GOBP_T_CELL_ACTIVATION.v2023.1.Hs.tsv")
pathway.list <- list(T_act = T_act$Gene,IEL = IEL.markers)
library(GSVA)

#1. ssGSEA
{
  Bulk.gsea <- gsva(TG_exp,pathway.list,method="ssgsea")
  TG_Samples$group <- desired_branch[rownames(TG_Samples)]
  df <- TG_Samples
  df$clusters <- NA
  df$clusters[df$SampleType %in% c("GTEx_Transverse","TCGA_colon","GTEx_Sigmoid")] <- "Colon"
  df$clusters[df$SampleType %in% c("GTEx_SI")] <- "GTEx_SI"
  df$clusters[df$group==2 & df$SampleType=="TCGA_CA"] <- "CA_C2"
  df$clusters[df$group==1 & df$SampleType=="TCGA_CA"] <- "CA_C1"
  df$clusters[df$group==1 & df$clusters=="Colon"] <- "Colon_C1"
  df$clusters[df$group==2 & df$clusters=="Colon"] <- "Colon_C2"
  df$IEL_Score <- scale(Bulk.gsea["IEL",])[,1]
  df$T_act <- scale(Bulk.gsea["T_act",])[,1]
  
  p1 <- ggplot(df,aes(clusters, IEL_Score)) + #stat_compare_means()+
    geom_boxplot(aes(fill = clusters),outlier.shape = NA,alpha = 1)+
    theme_classic()+ylab("IEL score")
  p2 <- ggplot(df,aes(clusters, T_act)) + #stat_compare_means()+
    geom_boxplot(aes(fill = clusters),outlier.shape = NA,alpha = 1)+
    theme_classic()+ylab("T_act score")
}

#2. CD160 expression across groups
{
  df$CD160 <- TG_exp["CD160",rownames(df)]
  p2 <- ggplot(df,aes(clusters, CD160)) + #stat_compare_means()+
    geom_boxplot(aes(fill = clusters),outlier.shape = NA,alpha = 1)+
    theme_classic()+ylab("CD160 Exp")
}
#----------------------------------------------------------------------------------------------------------
#### End ####
#----------------------------------------------------------------------------------------------------------













