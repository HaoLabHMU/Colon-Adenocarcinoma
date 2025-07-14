library(Seurat);library(dplyr);library(ggplot2);library(pheatmap);library(harmony)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer)
#----------------------------------------------------------------------------------------------------------
#### load data and functions ####
#----------------------------------------------------------------------------------------------------------
object <-  readRDS("allCells.rds")
load("scores.rda")
colour_bk <- c("#f0f0f0",colorRampPalette(c("#313695","#abd9e9"))(10),
               colorRampPalette(c("#abd9e9","#fee090"))(15),
               colorRampPalette(c("#fee090","#a50026"))(25))
ColAssign <- function(Var,palettes="Classic 20", n = 20){
  require(ggthemes);require(RColorBrewer)
  pal <- tableau_color_pal(palette = palettes,direction = 1,type="regular")
  if (length(Var) > n) {
    palOut <- colorRampPalette(pal(n))(length(Var))
    names(palOut) <- Var
  } else if (length(Var) == n) {
    palOut <- pal(n)
    names(palOut) <- Var
  } else if (length(Var) < n) {
    palOut <- pal(n)
    palOut <- setdiff(palOut,c("#7f7f7f","#c7c7c7"))# remove grey colors
    #palOut <- sample(palOut)
    palOut <- c(palOut,c("#7f7f7f","#c7c7c7"))
    palOut <- palOut[1:length(Var)]
    names(palOut) <- Var
  }
  return(palOut)
}
scCluster <- function(obj,nfeature=2500,min.d=0.1,res=1,TCRBCR=TRUE) {
  if (TCRBCR) {
    TCR.genes <- grep("^TR[AB][VJ]",rownames(obj),value = T)# keep GD TCR genes for GDT cells
    BCR.genes <- c(grep("^IG[KHL][VJC]",rownames(obj),value = T),
                   grep("^AC[0-9]",rownames(obj),value = T))# some RNA genes are also excluded.
  } else {TCR.genes <- c();BCR.genes <- c()}
  obj@assays$RNA@scale.data <- matrix() # in order to subset memory efficiently
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj,nfeatures = nfeature)
  var.genes <- VariableFeatures(obj)
  var.genes <- setdiff(var.genes,c(TCR.genes,BCR.genes))
  VariableFeatures(obj) <- var.genes
  obj <- ScaleData(obj,features = var.genes)
  obj <- RunPCA(obj, verbose = FALSE,features = VariableFeatures(obj),npcs=30)
  obj <- RunHarmony(obj, group.by.vars=c("orig.ident","Tissue2"),assay.use ="RNA")
  obj <- FindNeighbors(obj, dims=1:30,reduction = "harmony",k.param = 30)
  obj <- FindClusters(obj,resolution=res,random.seed=123,graph.name = 'RNA_snn')
  obj <- RunUMAP(obj,reduction = "harmony",seed.use = 123,dims=1:30,
                 umap.method='uwot',min.dist=min.d,spread=1)
  return(obj)
}
load("./Seurat_objects_afterN/color for CellType_n.rda")
load("./Seurat_objects_afterN/color for CellType_12.rda")
Tissue.colors <- c(PBMC="#F47D2B",LN="#f0f0f0",I="#33a02c",Colon="#89288F",T="#8A9FD1")
groupMeans <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
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

#----------------------------------------------------------------------------------------------------------
#### Figure 3A and related figures ####
#----------------------------------------------------------------------------------------------------------
metaData_fil <- object@meta.data
metaData_fil <- cbind(metaData_fil,scores[,-1])
# Natural.IEL
{
  IEL <- subset(CD8T_Tissue,subset=CellType_n=="Natural.IEL")
  IEL$CD160 <- GetAssayData(IEL,slot = "data")["CD160",]
  IEL$CD160 <- ifelse(IEL$CD160>0,"CD160+","CD160-")
  IEL$Cyt <- metaData_fil[colnames(IEL),"Cyt.score"]
  p1 <- ggplot(IEL@meta.data,aes(CD160, Cyt)) + #stat_compare_means()+
    geom_boxplot(aes(fill = Cyt),outlier.shape = NA,alpha = 1, col = "grey")+
    theme_classic()+ylab("Cytolytic activity")
}

# Act-T_ISG
{
  obj <- subset(CD8T_Tissue,subset=CellType_n=="CD8act_IFI")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(aes(fill = Exh.score),outlier.shape = NA,alpha = 1, col = "grey")+
    theme_classic()+ylab("Exh.score")
}

# Induced.IEL
{
  obj <- subset(CD8T_Tissue,subset=CellType_n=="Induced.IEL")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["Induced.IEL"])+
    theme_classic()+ylab("Exh.score")
}

# GZMK+ effector
{
  obj <- subset(CD8T_Tissue,subset=CellType_n=="GZMK+ effector")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["GZMK+ effector"])+
    theme_classic()+ylab("Exh.score")
}

# CD8_Mem
{
  obj <- subset(CD8T_Tissue,subset=CellType_n=="CD8_Mem")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["CD8_Mem"])+
    theme_classic()+ylab("Exh.score")
}

#----------------------------------------------------------------------------------------------------------
#### End ####
#----------------------------------------------------------------------------------------------------------
