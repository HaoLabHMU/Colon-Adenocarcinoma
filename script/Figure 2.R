library(Seurat);library(dplyr);library(ggplot2)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer)
#----------------------------------------------------------------------------------------------------------
#### load data and functions ####
#----------------------------------------------------------------------------------------------------------
object <-  readRDS("allCells.rds")
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
Tissue.colors <- c(PBMC="#F47D2B",LN="#f0f0f0",I="#33a02c",Colon="#89288F",T="#8A9FD1")
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
#----------------------------------------------------------------------------------------------------------
#### Fig.1B ####
#----------------------------------------------------------------------------------------------------------
cellType_1.colors <- ColAssign(unique(object$cellType_1),palettes="Tableau 10",n=10)
p1 <- DimPlot(object, reduction = "umap",label =T,group.by = "cellType_1",
              pt.size = 1,raster=T,shuffle=T,cols=cellType_1.colors)
pdf("./Analysis_afterN/Figures/All/umap by cellType_1.pdf",width = 6,height=6)
print(p1)
dev.off()
#----------------------------------------------------------------------------------------------------------
#### Fig.1C ####
#----------------------------------------------------------------------------------------------------------
library('ggsignif');library(viridis)
metaData <- object@meta.data
umaps <- Embeddings(object,reduction = "umap")
metaData <- cbind(metaData,umaps)
metaData_H <- metaData[metaData$Source=="Healthy",]
metaData_CA <- metaData[metaData$Source=="CA",]
for (i in unique(metaData_H$Tissue2)) {
  x <- metaData_H[metaData_H$Tissue2 == i,]
  g.overlay <- ggplot(data = x,aes(x = UMAP_1, y = UMAP_2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(color = 'white',size = .005)+
    scale_fill_viridis(option="D")+theme_classic()
  plot.path=paste0("./Analysis_afterN/Figures/All/umap by Density of ",i,"_healthy.pdf")
  ggsave(plot.path,width = 6,g.overlay)
}

for (i in unique(metaData_CA$Tissue2)) {
  x <- metaData_CA[metaData_CA$Tissue2 == i,]
  g.overlay <- ggplot(data = x,aes(x = UMAP_1, y = UMAP_2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
    geom_point(data = x[sample(seq_along(x[,1]),10000,replace = F),],
               aes(x = UMAP_1, y = UMAP_2),color = 'white',size = .005)+
    scale_fill_viridis(option="A")+theme_classic()
  plot.path=paste0("./Analysis_afterN/Figures/All/umap by Density of ",i,"_CA.pdf")
  ggsave(plot.path,width = 6,g.overlay)
}
#----------------------------------------------------------------------------------------------------------
#### Fig.1D ####
#----------------------------------------------------------------------------------------------------------
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)
library(ggalluvial)
load("color for CellType_12.rda")
df <- object@meta.data
df[,c("cellType_2")] <- factor(df[,c("cellType_2")])
df <- df %>% group_by(Source,Tissue2,cellType_2,.drop = FALSE) %>% summarise(n=n()) %>%mutate(freq = n/sum(n))
df$Tissue2 <- factor(df$Tissue2,levels=c("PBMC","LN","I","Colon","T"))
df$cellType_2 <- factor(df$cellType_2,levels=c('B cells','Plasma','CD4+ cells','CD8+ cells',
                                               'IEL','Other T','NK','MAST','Myeloid','pDC'))
# iOrd <- df$cells[df$CellType_n %in% c("Induced.IEL","Natural.IEL") & df$cellType_2=="Other T"]
# These 727 IEL-Ts are proliferating IEL-Ts, and all the proliferating cells are in "Other T" group.
p1 <- ggplot(df,aes(x = Tissue2,y = freq,alluvium =cellType_2,fill=cellType_2,
                    stratum = cellType_2)) +
  geom_alluvium(curve_type="linear",alpha = 0.9,knot.pos = 0.01) +
  geom_stratum(alpha = 1)+
  scale_x_discrete() +
  theme_classic()+
  scale_fill_manual("legend", values = cellType_2.colors) +
  facet_wrap(~Source) +
  ggtitle("alluvial plot for cellType_2")
#----------------------------------------------------------------------------------------------------------
#### Fig.1E ####
#----------------------------------------------------------------------------------------------------------
df <- object@meta.data
only_draw_pair_boxplot <- function(plot_df){
  p1 <- ggplot(data=plot_df, aes(x = Tissue2, y = value)) +
    geom_boxplot(alpha =0.7,size=1,outlier.shape = NA, mapping = aes(fill = Tissue2))+
    geom_jitter(size=3, shape=16,aes(group=Patient,col = Tissue2),
                alpha = 0.9,position = position_dodge(0))+
    scale_color_manual(values = Tissue.colors)+
    scale_fill_manual(values = Tissue.colors)+
    geom_line(aes(group = Patient), color = 'grey40', lwd = 0.3,position = position_dodge(0))+ #添加连线
    facet_wrap(~Source,ncol = 2)+
    theme_classic()
  return(p1)
}
detach("package:plyr", unload=TRUE); detach("package:reshape2", unload=TRUE)
df <- df %>% group_by(Source,Patient,Tissue2,cellType_1,.drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
df<-df[df$cellType_1=="Plasma",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
df$value <- df$freq
p1 <- only_draw_pair_boxplot(df)#fraction of Plasma cells

df <- object@meta.data
df$CD160 <- GetAssayData(object)["CD160",]
df$CD160 <- ifelse(df$CD160>0,"CD160+","CD160-")
df <- df %>% group_by(Source,Patient,Tissue2,CD160,.drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
df<-df[df$CD160=="CD160+",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
df$value <- df$freq
p1 <- only_draw_pair_boxplot(df)#fraction of CD160+ cells

df <- object@meta.data
df$NKG7 <- GetAssayData(object)["NKG7",]
df$NKG7 <- ifelse(df$NKG7>0,"NKG7+","NKG7-")
df <- df %>% group_by(Source,Patient,Tissue2,NKG7,.drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
df<-df[df$NKG7=="NKG7+",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
df$value <- df$freq
p1 <- only_draw_pair_boxplot(df)#fraction of NKG7+ cells
#----------------------------------------------------------------------------------------------------------
#### Fig.1J ####
#----------------------------------------------------------------------------------------------------------
library("factoextra");library(ggfortify);library(dendextend)
load("./Seurat_objects_afterN/color for CellType_12.rda")
# hierarchical clustering:
for (i in c("T cells","Plasma","NK","Myeloid","B cells","MAST")) {
  obj <- subset(object,subset=cellType_1==i)
  iOrd <- table(obj$orig.ident)
  iOrd <- names(iOrd)[iOrd>=50]
  obj <- subset(obj,subset=orig.ident %in% iOrd)
  obj <- FindVariableFeatures(obj,nfeatures = 2000)
  expMat <- GetAssayData(obj,slot = "data")[VariableFeatures(obj),]
  expMat <- groupMeans(expMat,groups=obj$orig.ident,sparse = T)
  
  dist_mat <- dist(t(expMat), method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'ward.D')
  plot(hclust_avg)
  avg_dend_obj <- as.dendrogram(hclust_avg)
  avg_col_dend <- color_branches(avg_dend_obj, k = 3)
  plot.path=paste0("hcluster by ",i,".pdf")
  pdf(plot.path,height = 4,width = 6)
  plot(avg_col_dend)
  dev.off()
}
#----------------------------------------------------------------------------------------------------------
#### Fig.1K ####
#----------------------------------------------------------------------------------------------------------
load("scores.rda")
df<- object@meta.data
df <- cbind(df,scores[,-1])
df<- df[df$cellType_1=="T cells",]
df$Exh.score
df$Tissue2 <- factor(df$Tissue2,levels=c("PBMC","LN","I","Colon","T"))

# Exhaustion score
ggplot(df,aes(Tissue2, Exh.score)) + 
  geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
  scale_fill_manual(values=Tissue.colors)+
  theme_classic()

# Exhaustion score
ggplot(df,aes(Tissue2, IEL.score)) + 
  geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
  scale_fill_manual(values=Tissue.colors)+
  theme_classic()

# Cytolytic score of IELs
df<- df[df$cellType_2=="IEL",]
ggplot(df,aes(Tissue2, Cyt.score)) + 
  geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
  scale_fill_manual(values=Tissue.colors)+
  theme_classic()
#----------------------------------------------------------------------------------------------------------
#### Fig.1M ####
#----------------------------------------------------------------------------------------------------------
# remove ribosomal genes, MT genes and BCR variable genes from analysis:
ribosomals <- grep(pattern = "^RP[SL]", x = rownames(object), value = TRUE)
TCR.genes <- grep("^TR[AB][VJ]",rownames(object),value = T)# keep GD TCR genes for GDT cells
BCR.genes <- c(grep("^IG[KHL][VJC]",rownames(object),value = T))
MT.genes <- grep(pattern = "^MT-", x = rownames(object), value = TRUE)

obj <- subset(object,subset = Source == "Healthy" & cellType_2 == "IEL")
Idents(obj) <- obj$Tissue2
DEGs <- FindMarkers(obj,ident.1 = "I",ident.2 = "Colon",logfc.threshold=0.1)
DEGs1 <- DEGs[rownames(DEGs) %ni% c(ribosomals,TCR.genes,MT.genes,BCR.genes),]

obj <- subset(object,subset = Source == "CA" & cellType_2 == "IEL")
Idents(obj) <- obj$Tissue2
DEGs <- FindMarkers(obj,ident.1 = "I",ident.2 = "Colon",logfc.threshold=0.1)
DEGs2 <- DEGs[rownames(DEGs) %ni% c(ribosomals,TCR.genes,MT.genes,BCR.genes),]

colnames(DEGs1) <- paste0(colnames(DEGs1),"-","Healthy")
colnames(DEGs2) <- paste0(colnames(DEGs2),"-","CA")
iOrd <- intersect(rownames(DEGs1),rownames(DEGs2))
df <- cbind(DEGs1[iOrd,],DEGs2[iOrd,]) # keep common DEGs
df <- df[rownames(df) != "MT2A",]#remove outliers
df$label[abs(df$`avg_log2FC-Healthy`)<0.5 | abs(df$`avg_log2FC-Healthy`)<0.5] <- NA

library(ggrepel)
ggplot(data=df, aes(x = `avg_log2FC-Healthy`, y = `avg_log2FC-CA`)) +
  geom_point(stroke = 0, alpha = 0.8,shape = 16,size=3) + 
  theme_classic() +
  scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")
#----------------------------------------------------------------------------------------------------------
#### End ####
#----------------------------------------------------------------------------------------------------------
