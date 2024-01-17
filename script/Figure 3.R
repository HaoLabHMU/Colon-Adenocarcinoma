library(Seurat);library(readxl);library(dplyr);library(ggplot2);library(pheatmap);library(harmony)
library(data.table);library(ggpubr);library(ArchR);library(DoubletFinder);library(RColorBrewer)
library(alakazam);
library(patchwork);detach("package:plyr", unload=TRUE); detach("package:reshape2", unload=TRUE)
setwd("~/Dropbox/Work/Harbin/ColonCancer")
load("./Seurat_objects_afterN/metaData_filtered.rda")
#----------------------------------------------------------------------------------------------------------
#### load data and functions ####
#----------------------------------------------------------------------------------------------------------
object <-  readRDS("allCells.rds")
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
B.colors <- c(B_Naive="#E7F6D2",B_Mem="#3EB7CC",Bfoc_MKI67="#ffd94a",
              Bfoc_NEIL1="#5F8E95",`Long-live PC`="#32A251",`Short-live PC`="#CE672E")
#----------------------------------------------------------------------------------------------------------
#### Figure 3A and S2A ####
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
#### Figure 3E, 3G and S2H ####
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
FeaturePlot(Bcell, c("IFITM2","IFITM1","STAT1","TGFB1"),raster=T,order=F,pt.size = 2) & 
  scale_colour_gradientn(colours = colour_bk)
#----------------------------------------------------------------------------------------------------------
#### Figure 3F, S2D and S2E ####
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
#### Figure 3H ####
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
#### Figure 3I ####
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
#### Figure S2J ####
#----------------------------------------------------------------------------------------------------------
load("monocle3 of PCs.rda")
df2 <- clone_mutation[clone_mutation$sequence_id %in% Bcell$cells,]
df2$cellType <- Bcell@meta.data[df2$sequence_id,"CellType_n"]
df2 <- df2[df2$cellType=="PCs",]
df2$pseudotime <- df[df2$sequence_id,'pseudotime']
p <- ggplot(df2, mapping = aes(x=pseudotime, y=100*mu_freq_seq_r)) + 
  theme_classic() + 
  xlab('Pseudotime') 
+ylab('SHM')
p <- p + 
  geom_smooth(aes(color = pseudotime), method = 'gam', se=T, color = 'black')
p1 <- p + 
  scale_y_continuous(breaks=seq(4, 6, 0.2))
#----------------------------------------------------------------------------------------------------------
#### End ####
#----------------------------------------------------------------------------------------------------------