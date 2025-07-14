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
B.colors <- c(B_Naive="#E7F6D2",B_Mem="#3EB7CC",Bfoc_MKI67="#ffd94a",
              Bfoc_NEIL1="#5F8E95",`Long-live PC`="#32A251",`Short-live PC`="#CE672E")
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
#### Extended Data Fig. 1A,1B ####
#----------------------------------------------------------------------------------------------------------
Bcell <- subset(object,subset= cellType_1 %in% c("B cells","Plasma"))

#1.clustering of B cells
{
  TCR.genes <- grep("^TR[AB][VJ]",rownames(Bcell),value = T)# keep GD TCR genes for GDT cells
  BCR.genes <- c(grep("^IG[KHL][VJC]",rownames(Bcell),value = T),
                 grep("^AC[0-9]",rownames(Bcell),value = T))# some RNA genes are also excluded.
  Bcell@assays$RNA@scale.data <- matrix() # in order to subset memory efficiently
  Bcell <- NormalizeData(Bcell)
  Bcell <- FindVariableFeatures(Bcell,nfeatures = 2500)
  var.genes <- VariableFeatures(Bcell)
  var.genes <- setdiff(var.genes,c(TCR.genes,BCR.genes))
  VariableFeatures(Bcell) <- var.genes
  Bcell <- ScaleData(Bcell,features = var.genes)
  Bcell <- RunPCA(Bcell, verbose = FALSE,features = VariableFeatures(Bcell),npcs=30)
  Bcell <- RunHarmony(Bcell, group.by.vars=c("orig.ident","Tissue2"),assay.use ="RNA")
  Bcell <- FindNeighbors(Bcell, dims=1:30,reduction = "harmony",k.param = 30)
  Bcell <- FindClusters(Bcell,resolution=1,random.seed=123,graph.name = 'RNA_snn')
  Bcell <- RunUMAP(Bcell,reduction = "harmony",seed.use = 123,dims=1:30,n.neighbors = 50,
                   umap.method='uwot',min.dist=0.8,spread=1)
}

#2. generate Umap and dotplot:
{
  load("monocle3 of PCs.rda")
  x <- df$cells[df$CellType_n=="Long-live"]
  y <- df$cells[df$CellType_n=="Short-live"]
  Bcell$cellType_2 <- paste0(Bcell$CellType_n)
  Bcell@meta.data[x,"cellType_2"] <- "Long-live PC"
  Bcell@meta.data[y,"cellType_2"] <- "Short-live PC"
  
  #1. UMAP
  B.colors <- c(B_Naive="#E7F6D2",B_Mem="#3EB7CC",Bfoc_MKI67="#ffd94a",
                Bfoc_NEIL1="#5F8E95",`Long-live PC`="#32A251",`Short-live PC`="#CE672E")
  p1 <- DimPlot(Bcell, reduction = "umap",label =F,group.by = "cellType_2",
                pt.size = 1,raster=T,shuffle=T,cols=B.colors)
  
  #Dotplot:
  Markers <- c("CD79A","CD19",#common markers
               "MS4A1","BANK1", "CCR7",#B_Mem
               "FCER2","TCL1A","IL4R","BACH2","IGHD","SELL",#B_Naive
               "NEIL1","RGS13","MEF2B","BCL6",'TUBA1B','MKI67','UBE2C','AURKB',#B_Foc
               "MZB1", "CD27",'JCHAIN',"IGHA1","IGHA2","STAT3","IKZF3"
  )
  DotPlot(Bcell, features = Markers,assay = "RNA",scale = T,group.by = "cellType_2") +
    scale_colour_gradientn(colors=rev(brewer.pal(9, "RdBu"))) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(breaks=Markers,labels=Markers)+ ylab("cellType_2")+ xlab("")
}
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 2A ####
#----------------------------------------------------------------------------------------------------------
DEGs <- FindMarkers(Bcell,ident.1 = "T",logfc.threshold=0.1,group.by = "Tissue2",min.pct = 0.3)
DEGs <- DEGs[rownames(DEGs) %ni% c(ribosomals,TCR.genes,MT.genes,BCR.genes),]
DEGs$gene <- rownames(DEGs)
DEGs <- DEGs[order(DEGs$avg_log2FC),]
write.table(DEGs,file = "DEGs of Tumor vs. other of B cells.txt",
            quote = F,sep = "\t",row.names = F) # DEGs for enrichment analysis of functions.

colour_bk <- c("#f0f0f0",colorRampPalette(c("#006837","#d9ef8b"))(5),
               colorRampPalette(c("#d9ef8b","#fee08b"))(5),
               colorRampPalette(c("#fee08b","#a50026"))(15))
FeaturePlot(Bcell, "TGFB1",raster=T,order=F,pt.size = 2) & 
  scale_colour_gradientn(colours = colour_bk)
#----------------------------------------------------------------------------------------------------------
#### Fig. 1F, Extended Data Fig. 2B,2C ####
#----------------------------------------------------------------------------------------------------------
metaData <- Bcell@meta.data
df <- metaData %>% group_by(Tissue2,Source,cellType_2,.drop = FALSE) %>% summarise(n=n())
df$cellType_2 <-factor(df$cellType_2, levels = names(B.colors))
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
p <- ggplot(df, aes(x = Tissue2, y = n, fill = cellType_2)) + 
  geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
  scale_fill_manual("legend", values = B.colors) +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(cols = vars(Source))#6*12

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
metaData <- Bcell@meta.data
df <- metaData %>% group_by(Source,Patient,Tissue2,cellType_2,.drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
df<-df[df$cellType_2=="Long-live PC",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
df$value <- df$freq
p1 <- only_draw_pair_boxplot(df)

load("clone_mutation.rda")
df <- clone_mutation %>% group_by(Source,Tissue2,Patient,isotype,.drop = FALSE) %>% summarise(n=n()) %>% mutate(freq = n/sum(n))
df<-df[df$isotype=="IgG",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
df$value <- df$freq
p2 <- only_draw_pair_boxplot(df)

df <- clone_mutation %>% group_by(Tissue2,Source,c_call,.drop = FALSE) %>% summarise(n=n())
df$c_call <-factor(df$c_call, levels = names(isotype.colors))
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
p3 <- ggplot(df, aes(x = Tissue2, y = n, fill = c_call)) + 
  geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
  scale_fill_manual("legend", values = isotype.colors) +
  theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(cols = vars(Source))#6*12
#----------------------------------------------------------------------------------------------------------
#### Figure 1G ####
#----------------------------------------------------------------------------------------------------------
#1. get clones from tumor B cells
load("clone_mutation.rda")
BCR.clone <- clone_mutation
df <- BCR.clone[BCR.clone$Tissue2=="T" ,]
metaData_fil <- Bcell@meta.data
metaData_fil <- metaData_fil[metaData_fil$cellType_1 %in% c("Plasma"),]
iOrd <- intersect(df$sequence_id,metaData_fil$cells)
df <- df[iOrd,]

#2. Get clones shared with normal tissues
df2 <- BCR.clone[BCR.clone$Tissue2 %in% c("I","Colon"),]
iOrd <- intersect(df2$sequence_id,metaData_fil$cells)
df2 <- df2[iOrd,]
iOrd <- intersect(df2$clone_id,df$clone_id)
df2 <- df2[df2$clone_id %in% iOrd,]
df <- df[df$clone_id %in% iOrd,]

#3. Pairwise clonal comparison:
df <- rbind(df,df2)
only_draw_pair_boxplot <- function(plot_df){
  p1 <- ggplot(data=plot_df, aes(x = Tissue2, y = value)) +
    geom_boxplot(alpha =0.7,size=1,outlier.shape = NA, mapping = aes(fill = Tissue2))+
    geom_jitter(size=3, shape=16,aes(group=clone_id,col = Tissue2),
                alpha = 0.9,position = position_dodge(0))+
    scale_color_manual(values = c(I="#33a02c",Colon="#89288F",T="#8A9FD1"))+
    scale_fill_manual(values = c(I="#33a02c",Colon="#89288F",T="#8A9FD1"))+
    geom_line(aes(group = clone_id), color = 'grey40', lwd = 0.3,position = position_dodge(0))+ #添加连线
    facet_wrap(~Source,ncol = 2)+
    theme_classic()
  return(p1)
}
df$TGFB1 <- GetAssayData(Bcell)["TGFB1",df$sequence_id]
df2 <- df %>% group_by(clone_id,Tissue2,Source,.drop = FALSE) %>% summarise(value=mean(TGFB1))
df2$Tissue2 <- factor(df2$Tissue2,levels=c("I","Colon","T"))
p1 <- only_draw_pair_boxplot(df2)
#----------------------------------------------------------------------------------------------------------
#### Figure 1H ####
#----------------------------------------------------------------------------------------------------------
df <- clone_mutation[clone_mutation$sequence_id %in% Bcell$cells,]
df$cellType <- Bcell@meta.data[df$sequence_id,"CellType_n"]
df <- df[df$cellType=="PCs",]
df <- df[df$Tissue2 %in% c("I","T","Colon"),]
df <-  df %>% group_by(patient,Tissue2,Source) %>% summarise(mut_freq_r=mean(mu_freq_seq_r))
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
df$value <- df$mut_freq_r
df$Patient <- df$patient
df$Tissue2 <- factor(df$Tissue2,levels = names(Tissue.colors))
p1 <- only_draw_pair_boxplot(df)

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
