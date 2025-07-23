library(Seurat);library(dplyr);library(ggplot2);library(pheatmap);library(harmony)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer)
#----------------------------------------------------------------------------------------------------------
#### load data and functions ####
#----------------------------------------------------------------------------------------------------------
object <-  readRDS("allCells.rds")
load("./tcr_data.rda")
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
  IEL <- subset(object,subset=CellType_n=="Natural.IEL")
  IEL$CD160 <- GetAssayData(IEL,slot = "data")["CD160",]
  IEL$CD160 <- ifelse(IEL$CD160>0,"CD160+","CD160-")
  IEL$Cyt <- metaData_fil[colnames(IEL),"Cyt.score"]
  p1 <- ggplot(IEL@meta.data,aes(CD160, Cyt)) + #stat_compare_means()+
    geom_boxplot(aes(fill = Cyt),outlier.shape = NA,alpha = 1, col = "grey")+
    theme_classic()+ylab("Cytolytic activity")
}

# Act-T_ISG
{
  obj <- subset(object,subset=CellType_n=="CD8act_IFI")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(aes(fill = Exh.score),outlier.shape = NA,alpha = 1, col = "grey")+
    theme_classic()+ylab("Exh.score")
}

# Induced.IEL
{
  obj <- subset(object,subset=CellType_n=="Induced.IEL")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["Induced.IEL"])+
    theme_classic()+ylab("Exh.score")
}

# GZMK+ effector
{
  obj <- subset(object,subset=CellType_n=="GZMK+ effector")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["GZMK+ effector"])+
    theme_classic()+ylab("Exh.score")
}

# CD8_Mem
{
  obj <- subset(object,subset=CellType_n=="CD8_Mem")
  obj$CD160 <- GetAssayData(obj,slot = "data")["CD160",]
  obj$CD160 <- ifelse(obj$CD160>0,"CD160+","CD160-")
  obj$Exh.score <- metaData_fil[colnames(obj),"Exh.score"]
  
  p1 <- ggplot(obj@meta.data,aes(CD160, Exh.score)) + #stat_compare_means()+
    geom_boxplot(outlier.shape = NA,alpha = 1,outlier.stroke=0,fill=T.colors["CD8_Mem"])+
    theme_classic()+ylab("Exh.score")
}
#----------------------------------------------------------------------------------------------------------
#### Figure 3B ####
#----------------------------------------------------------------------------------------------------------
df <- object@meta.data
df <- df[df$CellType_n=="Induced.IEL",]
iOrd <- intersect(df$cells,cell2clone.new$barcode)
df <- df[iOrd,];cell2clone.new <- cell2clone.new[iOrd,]
df$Patient_cloneID <- cell2clone.new[,'Patient_cloneID']
only_draw_pair_boxplot <- function(plot_df){
  p1 <- ggplot(data=plot_df, aes(x = Tissue2, y = value)) +
    geom_boxplot(alpha =0.7,size=1,outlier.shape = NA, mapping = aes(fill = Tissue2))+
    geom_jitter(size=3, shape=16,aes(group=Patient_cloneID,col = Tissue2),
                alpha = 0.9,position = position_dodge(0))+
    scale_color_manual(values = c(T="#8A9FD1",Colon="#89288F",I="#389937"))+
    scale_fill_manual(values = c(T="#8A9FD1",Colon="#89288F",I="#389937"))+
    geom_line(aes(group = Patient_cloneID), color = 'grey40', lwd = 0.3,position = position_dodge(0))+ #添加连线
    facet_wrap(~Source,ncol = 2)+
    theme_classic()
  return(p1)
}

#keep clones across tissues
df2 <- df %>% group_by(Patient_cloneID,Tissue2) %>% summarise(n=n())
iOrd <- table(df2$Patient_cloneID)
iOrd <- names(iOrd)[iOrd>1]
df <- df[df$Patient_cloneID %in% iOrd,]
df$Tex.score <- scores[df$cells,"Exh.score"]
df2 <- df %>% group_by(Patient_cloneID,Tissue2,Source,.drop = FALSE) %>% summarise(Tex.score=mean(Tex.score))
df2$Tissue2 <- factor(df2$Tissue2,levels=c("I","Colon","T"))
df2$value <- df2$Tex.score
p1 <- only_draw_pair_boxplot(df2)


Prog.Exh <- c("IL7R","GPR183","LMNA","NR4A3","TCF7","MGAT4A",
              "CD55","AIM1","PER1","FOSL2","EGR1","TSPYL2",
              "YPEL5","CSRNP1","REL","SKIL","PIK3R1","FOXP1",
              "RGCC","PFKFB3","MYADM","ZFP36L2","USP36",
              "TC2N","FAM177A1","BTG2","TSC22D2","FAM65B","STAT4",
              "RGPD5","NEU1","IFRD1","PDE4B","NR4A1")
Tcell <- subset(object, subset=cellType_1=="T cells")
cells_AUC <- AUCell_run(Tcell@assays$RNA@counts, Prog.Exh,aucMaxRank=3000)
Prog.scores <- getAUC(cells_AUC)["geneSet",]
df$Pex.score <- Prog.scores[df$cells]
df2 <- df %>% group_by(Patient_cloneID,Tissue2,Source,.drop = FALSE) %>% summarise(Pex.score=mean(Pex.score))
df2$Tissue2 <- factor(df2$Tissue2,levels=c("I","Colon","T"))
df2$value <- df2$Pex.score
p1 <- only_draw_pair_boxplot(df2)

#----------------------------------------------------------------------------------------------------------
#### Figure 3C ####
#----------------------------------------------------------------------------------------------------------
load("./tcr_data.rda")
#1. get clones from Induced.IEL from tumor
df <- object@meta.data
df <- df[df$CellType_n=="Induced.IEL" & df$Tissue2=="T",]
iOrd <- intersect(df$cells,cell2clone.new$barcode)
df <- df[iOrd,];Clones <- cell2clone.new[iOrd,]
df$Patient_cloneID <- Clones[,'Patient_cloneID']

#2. Get clones shared with LN or PBMC
df.PLN <- object@meta.data
df.PLN <- df.PLN[df.PLN$Tissue2 %in% c("LN","PBMC"),]
iOrd <- intersect(df.PLN$cells,cell2clone.new$barcode)
df.PLN <- df.PLN[iOrd,];Clones <- cell2clone.new[iOrd,]
df.PLN$Patient_cloneID <- Clones[,'Patient_cloneID']
iOrd <- intersect(df.PLN$Patient_cloneID,df$Patient_cloneID)
df.PLN <- df.PLN[df.PLN$Patient_cloneID %in% iOrd,]
df <- df[df$Patient_cloneID %in% iOrd,]

#3. Pairwise clonal comparision:
{
  iOrd <- intersect(colnames(df),colnames(df.PLN))
  df <- rbind(df[,iOrd],df.PLN[,iOrd])
  only_draw_pair_boxplot <- function(plot_df){
    p1 <- ggplot(data=plot_df, aes(x = Tissue2, y = value)) +
      geom_boxplot(alpha =0.7,size=1,outlier.shape = NA, mapping = aes(fill = Tissue2))+
      geom_jitter(size=3, shape=16,aes(group=Patient_cloneID,col = Tissue2),
                  alpha = 0.9,position = position_dodge(0))+
      scale_color_manual(values = c(PBMC="#F47D2B",LN="#f0f0f0",T="#8A9FD1"))+
      scale_fill_manual(values = c(PBMC="#F47D2B",LN="#f0f0f0",T="#8A9FD1"))+
      geom_line(aes(group = Patient_cloneID), color = 'grey40', lwd = 0.3,position = position_dodge(0))+ #添加连线
      facet_wrap(~Source,ncol = 2)+
      theme_classic()
    return(p1)
  }
  
  # Terminal Exhaustion score:
  {
    df$Tex.score <- scores[df$cells,"Exh.score"]
    df2 <- df %>% group_by(Patient_cloneID,Tissue2,Source,.drop = FALSE) %>% summarise(Tex.score=mean(Tex.score))
    df2$Tissue2 <- factor(df2$Tissue2,levels=c("PBMC","LN","T"))
    df2$value <- df2$Tex.score
    
    p1 <- only_draw_pair_boxplot(df2)
  }
  
  # Proginetor score:
  {
    df$Pex.score <- Prog.scores[df$cells]
    df2 <- df %>% group_by(Patient_cloneID,Tissue2,Source,.drop = FALSE) %>% 
      summarise(Pex.score=mean(Pex.score))
    df2$Tissue2 <- factor(df2$Tissue2,levels=c("PBMC","LN","T"))
    df2$value <- df2$Pex.score
    p1 <- only_draw_pair_boxplot(df2)
  }
}
#----------------------------------------------------------------------------------------------------------
#### Figure 3D ####
#----------------------------------------------------------------------------------------------------------
df <- object@meta.data
CD8T.cells <- df[df$cellType_2 %in% c("CD8+ cells","IEL","Other T"),]
CD8T.cells <- CD8T.cells[CD8T.cells$CellType_n != "Natural.IEL",]
# only keep cells from Tumor
CD8T.cells <- CD8T.cells[CD8T.cells$Tissue2 %in% c("T"),]
CD8T.cells$CD160 <- GetAssayData(object)["CD160",CD8T.cells$cells]
CD8T.cells$CD160 <- ifelse(CD8T.cells$CD160>0,"Positive","Negative")
CD8T.cells[,c("CD160")] <- factor(CD8T.cells$CD160)
CD8T.cells$Pex.scores <- Prog.scores[CD8T.cells$cells]
CD8T.cells$Tex.scores <- scores[CD8T.cells$cells,"Exh.score"]

boxplot(CD8T.cells$Pex.scores ~ CD8T.cells$CD160)
boxplot(CD8T.cells$Tex.scores ~ CD8T.cells$CD160)
#----------------------------------------------------------------------------------------------------------
#### Figure 3F ####
#----------------------------------------------------------------------------------------------------------
load("./tcr_data.rda")
df <- object@meta.data
CD8T.cells <- df[df$cellType_2 %in% c("CD8+ cells","IEL","Other T"),]
CD8T.cells <- CD8T.cells[CD8T.cells$CellType_n != "Natural.IEL",]

# only keep Ileum,Colon, Tumor and PBMC
CD8T.cells <- CD8T.cells[CD8T.cells$Tissue2 %in% c("I","Colon","T","PBMC"),]
CD8T.cells$CD160 <- GetAssayData(object)["CD160",CD8T.cells$cells]
CD8T.cells$CD160 <- ifelse(CD8T.cells$CD160>0,"Positive","Negative")
CD8T.cells[,c("CD160")] <- factor(CD8T.cells$CD160)

#keep CD8T cells matched with TCR:
iOrd <- intersect(CD8T.cells$cells,cell2clone.new$barcode)
df <- CD8T.cells[iOrd,]
df$clone_size <- cell2clone.new[iOrd,'Freq.Patient']
df$Tissue <- paste0(df$Source,"-",df$Tissue2)
df$clone_id <- cell2clone.new[iOrd,'Patient_cloneID']

df$group <- case_when(
  df$clone_size == 1 ~ "1",
  df$clone_size >= 2 & df$clone_size <= 9 ~ "2",
  df$clone_size >= 10 & df$clone_size <= 30 ~ "3",
  df$clone_size >= 31 ~ "4")
x <- table(df$group,df$Tissue)
y <- table(df$group[df$CD160=="Positive"],df$Tissue[df$CD160=="Positive"])
x <- melt(x); y <- melt(y)
colnames(x) <- c("Group","Tissue","N.Cells")

# Data for pie chat
piechat_frac <- x
piechat_frac$N.CD160Cells <- y$value
piechat_frac$Frac <- piechat_frac$N.CD160Cells/piechat_frac$N.Cells
#----------------------------------------------------------------------------------------------------------
#### End ####
#----------------------------------------------------------------------------------------------------------
