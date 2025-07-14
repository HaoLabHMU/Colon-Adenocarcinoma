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
#### Figure 2A ####
#----------------------------------------------------------------------------------------------------------
CD8.Tcell <- subset(object,subset = CellType_n %in% 
                      c('CD8_Mem','CD8_Naive','CD8act_IFI',
                        'GZMK+ effector','Induced.IEL',
                        'MAIT','Natural.IEL'))
CD8T_Tissue <- subset(CD8.Tcell,subset=Tissue2 %in% c("I","Colon","T"))
CD8T_Tissue <- scCluster(CD8T_Tissue,nfeature=2500,min.d=0.1,res=1.5)
CD8T_Tissue <- RunUMAP(CD8T_Tissue,reduction = "harmony",seed.use = 123,n.neighbors = 50,
                       dims=1:30,min.dist=0.6,spread=1)

p1 <- DimPlot(CD8T_Tissue, reduction = "umap",label =T,group.by = "CellType_n",
              pt.size = 1.5,raster=T,shuffle=T,cols=T.colors)

CD8T_Tissue$IEL.score <- scores[colnames(CD8T_Tissue),"IEL.score"]
p2 <- FeaturePlot(CD8T_Tissue, features = "IEL.score",
                  raster=T,reduction = "umap",pt.size = 2) & 
  scale_colour_gradientn(colours = colour_bk)

p3 <- FeaturePlot(CD8T_Tissue, features = "CD6",
                  raster=T,reduction = "umap",pt.size = 2)
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 3G ####
#----------------------------------------------------------------------------------------------------------
sel.genes <- c("CD8A",
               "LEF1","TCF7","CCR7","SELL",#Naive
               "ANXA1","IL7R","LMNA","LTB", #Tcm
               "GZMK","GZMA","GZMB","PRF1","NKG7","GZMH", #Tem and CTL
               "IFIT1","IFIT2","IFIT3",# IFI
               "RORC","SLC4A10","TRAV1-2","CCR6",#MAIT
               IEL.markers,
               "TRDC","TRDV1","TYROBP"
)
CD8T_Tissue$CellType_n <- factor(CD8T_Tissue$CellType_n,levels=c(
  'CD8_Naive','CD8_Mem','GZMK+ effector','CD8act_IFI','MAIT','Induced.IEL','Natural.IEL'))

p1 <- DotPlot(CD8T_Tissue, features = sel.genes,assay = "RNA",scale = T,group.by = "CellType_n") +
  scale_colour_gradientn(colors=rev(brewer.pal(9, "RdBu"))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=sel.genes,labels=sel.genes)+ ylab("Tissue")+ xlab("")
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 3H,3I ####
#----------------------------------------------------------------------------------------------------------
CD4T <- subset(object, subset=cellType_2=="CD4+ cells")
CD4T <- scCluster(CD4T,nfeature=2500,min.d=0.1,res=1)
CD4T <- RunUMAP(CD4T,reduction = "harmony",seed.use = 123,n.neighbors = 50,
                dims=1:30,min.dist=0.6,spread=1)
p1 <- DimPlot(CD4T, reduction = "umap",label =T,group.by = "CellType_n",
              pt.size = 1,raster=T,shuffle=T,cols=T.colors)

sel.genes <- c("LEF1","TCF7","CCR7","SELL",#Naive
               "ANXA1","IL7R","LMNA","CD69","LTB", #Tcm
               "GIMAP4","GZMK","GZMA","PRF1","NKG7","GZMH","GNLY", #Tem and CTL
               "CXCL13","BCL6","CXCR5","IL21","TOX2","TOX","PDCD1",#Tfh
               "FOXP3","IL2RA","CTLA4"#"Treg"
)
CD4T$CellType_n <- factor(CD4T$CellType_n,levels=c(
  'CD4_Naive','CD4_Tcm','CD4_Tem','CD4_CTL','Tfh','Treg'))
p2 <- DotPlot(CD4T, features = sel.genes,assay = "RNA",scale = T,group.by = "CellType_n") +
  scale_colour_gradientn(colors=rev(brewer.pal(9, "RdBu"))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=sel.genes,labels=sel.genes)+ ylab("Tissue")+ xlab("")
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 3J,3K ####
#----------------------------------------------------------------------------------------------------------
CD8.Tcell <- subset(object,subset = CellType_n %in% 
                      c('CD8_Mem','CD8_Naive','CD8act_IFI',
                        'GZMK+ effector','Induced.IEL',
                        'MAIT','Natural.IEL'))
PLN <- subset(CD8.Tcell,subset=Tissue2 %in% c("PBMC","LN"))
PLN <- scCluster(PLN,nfeature=2500,min.d=0.1,res=1.5)
PLN <- RunHarmony(PLN, group.by.vars=c("orig.ident","Tissue2"),assay.use ="RNA")
PLN <- FindNeighbors(PLN, dims=1:30,reduction = "harmony",k.param = 30)
PLN <- FindClusters(PLN,resolution=1,random.seed=123,graph.name = 'RNA_snn')
PLN <- RunUMAP(PLN,reduction = "harmony",seed.use = 123,dims=1:30,n.neighbors = 100,
               umap.method='uwot',min.dist=0.2,spread=1)
p1 <- DimPlot(PLN, reduction = "umap",label =T,group.by = "CellType_n",
              pt.size = 1.3,raster=T,shuffle=T,cols=T.colors)

sel.genes <- c("LEF1","TCF7","CCR7","SELL",#Naive
               "KLF2","IL7R","LMNA","LTB", #Tcm
               "GZMK","GZMA","GZMB","PRF1","NKG7","GNLY", #Tem and CTL
               "LAG3","TOX","PDCD1","HAVCR2",#Exh
               "RORC","SLC4A10","TRAV1-2","CEBPD","CCR6" #MAIT
)
PLN$CellType_n <- factor(PLN$CellType_n,levels=c(
  'CD8_Naive','CD8_Mem','CD8_CTL','CD8_Exh','MAIT'))
p2 <- DotPlot(PLN, features = sel.genes,assay = "RNA",scale = T,group.by = "CellType_n") +
  scale_colour_gradientn(colors=rev(brewer.pal(9, "RdBu"))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=sel.genes,labels=sel.genes)+ ylab("Tissue")+ xlab("")
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 1H ####
#----------------------------------------------------------------------------------------------------------
library("scales")
metaData <- object@meta.data[c(colnames(CD4T),colnames(PLN),colnames(CD8T_Tissue)),]
metaData$Tissue <- paste0(metaData$Source,"_",metaData$Tissue2)
x <- table(metaData$Tissue)/nrow(metaData)
expected <- as.numeric(x);names(expected) <- names(x)
expected <- expected[c('CA_Colon','CA_I','CA_LN','CA_PBMC','CA_T','Healthy_Colon','Healthy_I')]
df <- NULL
cellTypes <- unique(metaData$CellType_n)
for (i in cellTypes) {
  x <- metaData[metaData$CellType_n==i,]
  temp <- c(CA_Colon=0,CA_I=0,CA_LN=0,CA_PBMC=0,CA_T=0,Healthy_Colon=0,Healthy_I=0)
  temp[names(table(x$Tissue))] <- table(x$Tissue)
  x <- temp/nrow(x)/expected
  x <- data.frame(cellType=i,Tissue=names(x),Roe=as.numeric(x))
  df <- rbind(df,x)
}

labs <- as.character(round(df$Roe,2))
df2 <- df;df2$Roe[df2$Roe>=2.3]=2.3#change the max value to 2
iOrd <- c('CD4_Naive','CD4_Tcm','CD4_Tem','CD4_CTL','Tfh','Treg',
          'CD8_Naive','CD8_Mem','GZMK+ effector','CD8act_IFI','CD8_CTL','CD8_Exh',
          'MAIT','Induced.IEL','Natural.IEL')
iOrd2 <- c('CA_PBMC','CA_LN','CA_I','CA_Colon','CA_T','Healthy_I','Healthy_Colon')
df2$cellType <- factor(df2$cellType,levels = iOrd)
df2$Tissue <- factor(df2$Tissue,levels = iOrd2)
p1 <- ggplot(data =  df2, aes(x = Tissue, y = cellType)) + 
  geom_tile(aes(fill = Roe)) +
  scale_fill_gradientn(colours=c("#4393c3","#92c5de","#f7f7f7","#f4a582","#d6604d","#b30000"),
                       values=rescale(c(0,0.5,1,1.5,2,2.3)),
                       guide="colorbar")
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 1I ####
#----------------------------------------------------------------------------------------------------------
detach("package:plyr", unload=TRUE)
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

metaData_fil <- object@meta.data
df <- metaData_fil[metaData_fil$cellType_2 == "CD4+ cells",]
df <- df %>% group_by(Source,Patient,Tissue2,CellType_n,.drop = FALSE) %>% 
  summarise(n = n()) %>% mutate(freq = n/sum(n))

df2<-df[df$CellType_n=="Treg",]
df2$Tissue2 <-factor(df2$Tissue2, levels = names(Tissue.colors))
df2$value <- df2$freq
p1 <- only_draw_pair_boxplot(df2)

df2<-df[df$CellType_n=="Tfh",]
df2$Tissue2 <-factor(df2$Tissue2, levels = names(Tissue.colors))
df2$value <- df2$freq
p2 <- only_draw_pair_boxplot(df2)
#----------------------------------------------------------------------------------------------------------
#### Extended Data Fig. 3N,3O ####
#----------------------------------------------------------------------------------------------------------
#1. across cellTypes of CD8 T and CD4 T
{
  df <- object@meta.data[colnames(CD8T_Tissue),]
  df <- cbind(df,scores[df$cells,-1])
  df$CellType_n <- factor(df$CellType_n,levels=c(
    'CD8_Naive','CD8_Mem','GZMK+ effector','CD8act_IFI','MAIT','Induced.IEL','Natural.IEL'))
  p1 <- ggplot(df,aes(CellType_n, NeoTCR_CD8)) + 
    geom_boxplot(aes(fill = CellType_n),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=all.cols)+
    theme_classic()+ylab("Neoantigen reactivity (CD8+)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  df <- object@meta.data[colnames(CD4T),]
  df <- cbind(df,scores[df$cells,-1])
  df$CellType_n <- factor(df$CellType_n,levels=c(
    'CD4_Naive','CD4_Tcm','CD4_Tem','CD4_CTL','Tfh','Treg'))
  p2 <- ggplot(df,aes(CellType_n, NeoTCR_CD4)) + 
    geom_boxplot(aes(fill = CellType_n),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=all.cols)+
    theme_classic()+ylab("Neoantigen reactivity (CD4+)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#2. across Tissues for CD8 T and CD4 T
{
  df <- object@meta.data
  df <- cbind(df,scores[df$cells,-1])
  df <- df[df$CellType_n=="Treg",]
  df$Tissue2 <- factor(df$Tissue2,levels = names(Tissue.colors))
  p1 <- ggplot(df,aes(Tissue2, NeoTCR_CD4)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Neoantigen reactivity (CD4+)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~Source)
  
  df <- object@meta.data
  df <- cbind(df,scores[df$cells,-1])
  df <- df[df$CellType_n=="Tfh",]
  df$Tissue2 <- factor(df$Tissue2,levels = names(Tissue.colors))
  p2 <- ggplot(df,aes(Tissue2, NeoTCR_CD4)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Neoantigen reactivity (CD4+)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~Source)
}

#3. across Tissues for IEL-Ts
{
  df <- object@meta.data
  df <- cbind(df,scores[df$cells,-1])
  df <- df[df$CellType_n=="Induced.IEL",]
  df$Tissue2 <- factor(df$Tissue2,levels = names(Tissue.colors))
  p3 <- ggplot(df,aes(Tissue2, NeoTCR_CD8)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Neoantigen reactivity (CD8+)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~Source)
  
  p4 <- ggplot(df,aes(Tissue2, Exh.score)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Exh.score")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~Source)
  
  p5 <- ggplot(df,aes(Tissue2, NeoTCR_CD8)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Neoantigen reactivity (CD8)")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~Source)
  
  p6 <- ggplot(df,aes(Tissue2, Virus_Specific)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Virus Specificity")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~Source)
  
  
  df <- object@meta.data
  df <- cbind(df,scores[df$cells,-1])
  df <- df[df$CellType_n=="Natural.IEL",]
  df$Tissue2 <- factor(df$Tissue2,levels = names(Tissue.colors))
  p7 <- ggplot(df,aes(Tissue2, Cyt.score)) + 
    geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1)+
    scale_fill_manual(values=Tissue.colors)+
    theme_classic()+ylab("Cyt.score")+
    theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_wrap(~Source)
}
#----------------------------------------------------------------------------------------------------------
#### Figure 2B and related figures ####
#----------------------------------------------------------------------------------------------------------
sel.genes1 <- c("IFNG","NKG7","CGAS","RGS1","CCR9","IL7R","CD69","DUSP2")
sel.genes2 <- c("LYST","CTLA4","TIGIT","LAYN","CD74","PDCD1","LAG3","HAVCR2")
bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
               colorRampPalette(c("#fddbc7","#b2182b"))(84))
mat <- GetAssayData(CD4T)[c(sel.genes1,sel.genes2),]
Mat <- groupMeans(mat,groups=CD4T$Tissue2,na.rm = TRUE, sparse = T)
Mat <- Mat[,c("I","Colon","T")]
pheatmap(Mat, cluster_rows=F, cluster_cols=F,
         scale="row",show_colnames=T,show_rownames=T,
         breaks = bk,color = colour_bk,angle_col=90)

mat <- GetAssayData(CD8T_Tissue)[c(sel.genes1,sel.genes2),]
Mat <- groupMeans(mat,groups=CD8T_Tissue$Tissue2,na.rm = TRUE, sparse = T)
pheatmap(Mat, cluster_rows=F, cluster_cols=F,
         scale="row",show_colnames=T,show_rownames=T,
         breaks = bk,color = colour_bk,angle_col=90)

only_draw_pair_boxplot <- function(plot_df){
  p1 <- ggplot(data=plot_df, aes(x = Tissue2, y = value)) +
    geom_boxplot(alpha =0.7,size=1,outlier.shape = NA, mapping = aes(fill = Tissue2))+
    geom_jitter(size=3, shape=16,aes(group=Patient,col = Tissue2),
                alpha = 0.9,position = position_dodge(0))+
    scale_color_manual(values = Tissue.colors)+
    scale_fill_manual(values = Tissue.colors)+
    geom_line(aes(group = Patient), color = 'grey40', lwd = 0.3,position = position_dodge(0))+ #添加连线
    theme_classic()
  return(p1)
}
# CD4 T
{
  CD4T$class <- paste0(CD4T$Tissue2,"-",CD4T$Patient)
  mat <- GetAssayData(CD4T)[c(sel.genes1,sel.genes2),]
  Mat <- groupMeans(mat,groups=CD4T$class,na.rm = TRUE, sparse = T)
  df <- unique(CD4T@meta.data[,c("class","Tissue2","Patient")])
  rownames(df) <- df$class
  df <- df[colnames(Mat),]
  df$sel1 <- colMeans(Mat[sel.genes1,])
  df$sel2 <- colMeans(Mat[sel.genes2,])
  df$value = df$sel1
  df$Tissue2 <- factor(df$Tissue2,levels = c("I","Colon","T"))
  p1 <- only_draw_pair_boxplot(df)
  df$value = df$sel2
  p2 <- only_draw_pair_boxplot(df)
}
# CD8 T
{
  CD8T_Tissue$class <- paste0(CD8T_Tissue$Tissue2,"-",CD8T_Tissue$Patient)
  mat <- GetAssayData(CD8T_Tissue)[c(sel.genes1,sel.genes2),]
  Mat <- groupMeans(mat,groups=CD8T_Tissue$class,na.rm = TRUE, sparse = T)
  df <- unique(CD8T_Tissue@meta.data[,c("class","Tissue2","Patient")])
  rownames(df) <- df$class
  df <- df[colnames(Mat),]
  df$sel1 <- colMeans(Mat[sel.genes1,])
  df$sel2 <- colMeans(Mat[sel.genes2,])

  df$Tissue2 <- factor(df$Tissue2,levels = c("I","Colon","T"))
  df$value = df$sel1
  p1 <- only_draw_pair_boxplot(df)
  df$value = df$sel2
  p2 <- only_draw_pair_boxplot(df)
}
#----------------------------------------------------------------------------------------------------------
#### Figure 2G and related figures ####
#----------------------------------------------------------------------------------------------------------
library(monocle)
IEL <- subset(object,subset=CellType_n=="Induced.IEL")
IEL <- FindVariableFeatures(IEL,nfeatures = 2500)
cell_metadata <- IEL@meta.data
expression_data <- GetAssayData(IEL,slot="counts")
gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                              stringsAsFactors = F,row.names = rownames(expression_data))
pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
IEL.cds <- newCellDataSet(cellData=expression_data,phenoData = pd,
                          featureData = fd,expressionFamily=negbinomial.size())
IEL.cds <- estimateSizeFactors(IEL.cds)
IEL.cds <- estimateDispersions(IEL.cds)
ordering_genes <- VariableFeatures(IEL)
IEL.cds <- setOrderingFilter(IEL.cds, ordering_genes)
IEL.cds <- reduceDimension(IEL.cds, max_components = 2,reduction_method = 'DDRTree',pseudo_expr=1)
IEL.cds <- orderCells(IEL.cds)
IEL.cds$Exh.score <- scores[colnames(IEL.cds),"Exh.score"]
IEL.cds$IEL.score <- scores[colnames(IEL.cds),"IEL.score"]

df=phenoData(IEL.cds)@data
df <- cbind(df,t(IEL.cds@reducedDimS))#c("Component_1","Component_2")
df <- cbind(df,scores[rownames(df),c("NeoTCR_CD8","Tumor_Specific","Virus_Specific","NeoTCR_CD4","MANA_TIL","Influenza_TIL")])
iOrd <- c('MKI67','PDCD1',"CD160")
df <- cbind(df,t(GetAssayData(IEL,slot="data")[iOrd,rownames(df)]))
colour_bk <- c(colorRampPalette(c("#313695","#abd9e9"))(10),
               colorRampPalette(c("#abd9e9","#fee090"))(15),
               colorRampPalette(c("#fee090","#a50026"))(15))
p <- plot_cell_trajectory(IEL.cds, color_by = "Pseudotime",cell_size=1,show_branch_points=F)+
  scale_color_gradientn(colors=colour_bk)

colour_bk <- c("#f0f0f0",colorRampPalette(c("#313695","#abd9e9"))(10),
               colorRampPalette(c("#abd9e9","#fee090"))(15),
               colorRampPalette(c("#fee090","#a50026"))(15))
p1 <- ggplot(df, aes(x=Component_1, y=Component_2)) +
  geom_point(aes(color=MKI67),size=0.5)+
  scale_color_gradientn(colors=colour_bk)+
  theme_classic()
p2 <- ggplot(df, aes(x=Component_1, y=Component_2)) +
  geom_point(aes(color=Exh.score),size=0.5)+
  scale_color_gradientn(colors=colour_bk)+
  theme_classic()

p1 <- ggplot(df, mapping = aes(x=Pseudotime, y=Exh.score)) + 
  theme_classic() + xlab('Pseudotime') + ylab('Exh.score')+
  geom_smooth(aes(group = Tissue2,color = Tissue2),method = 'loess', se=F)+
  scale_color_manual(values=Tissue.colors)
p2 <- ggplot(df, mapping = aes(x=Pseudotime, y=NeoTCR_CD8)) + 
  theme_classic() + xlab('Pseudotime') +
  geom_smooth(aes(group = Tissue2,color = Tissue2),method = 'loess', se=F)+
  scale_color_manual(values=Tissue.colors)
p3 <- ggplot(df, mapping = aes(x=Pseudotime, y=Virus_Specific)) + 
  theme_classic() + xlab('Pseudotime') +
  geom_smooth(aes(group = Tissue2,color = Tissue2),method = 'loess', se=F)+
  scale_color_manual(values=Tissue.colors)
p4 <- ggplot(df, mapping = aes(x=Pseudotime, y=Tumor_Specific)) + 
  theme_classic() + xlab('Pseudotime') +
  geom_smooth(aes(group = Tissue2,color = Tissue2),method = 'loess', se=F)+
  scale_color_manual(values=Tissue.colors)

# gene expression across pseudotime
{
  for (i in genes){
    i=ensym(i)
    p <- ggplot(df, mapping = aes(x=Pseudotime, y = !!i)) + 
      theme_classic() + xlab('Pseudotime') + ylab('expression')+
      geom_smooth(aes(group = Tissue2,color = Tissue2),method = 'loess', se=F)+
      scale_color_manual(values=Tissue.colors)
    plot.path <- paste0("Smooth line of ",i," in IEL trajectory.pdf")
    pdf(plot.path,width = 6,height = 4)
    print(p)
    dev.off()
  }
}
#----------------------------------------------------------------------------------------------------------
#### Figure 2I ####
#----------------------------------------------------------------------------------------------------------
df <- object@meta.data[colnames(CD8.Tcell),]
df$CD160 <- GetAssayData(CD8.Tcell,slot="data")["CD160",]
df <- df[df$CellType_n != "Natural.IEL",]#remove GDT
df$CD160 <- ifelse(df$CD160>0,"Positive","Negative")
df[,c("CD160")] <- factor(df$CD160)

detach("package:plyr", unload=TRUE)
df <- df %>% group_by(Source,Patient,Tissue2,CD160,.drop = FALSE) %>% summarise(n = n()) %>% mutate(freq = n/sum(n))
df<-df[df$CD160=="Positive",]
df$Tissue2 <-factor(df$Tissue2, levels = names(Tissue.colors))
p1 <- ggplot(df,aes(Tissue2, freq)) + #stat_compare_means()+
  geom_boxplot(aes(fill = Tissue2),outlier.shape = NA,alpha = 1, col = "grey")+
  geom_jitter(width = 0.1)+scale_fill_manual(values=Tissue.colors)+
  theme_classic()+ylab("Fraction of CD160+ T cells")+
  facet_wrap(~Source)

#----------------------------------------------------------------------------------------------------------
#### Figure 2M and related figures ####
#----------------------------------------------------------------------------------------------------------
CD8T_Tissue$CD160 <- GetAssayData(CD8T_Tissue)["CD160",]
CD8T_Tissue$CD160 <- ifelse(CD8T_Tissue$CD160>0,"CD160+","CD160-")
CD8T_Colon <- subset(CD8T_Tissue, subset=Tissue2=="Colon")
DEGs_Colon <- FindMarkers(CD8T_Colon,ident.1 = "CD160+",ident.2="CD160-",group.by = "CD160",logfc.threshold = 0.1)
DEGs_Colon$genes <- rownames(DEGs_Colon)
CD8T_Ileum <- subset(CD8T_Tissue, subset=Tissue2=="I")
DEGs_Ileum <- FindMarkers(CD8T_Ileum,ident.1 = "CD160+",ident.2="CD160-",group.by = "CD160",logfc.threshold = 0.1)
DEGs_Ileum$genes <- rownames(DEGs_Ileum)

iOrd <- intersect(rownames(DEGs_Ileum),rownames(DEGs_Colon))# keep common DEGs
df <- DEGs_Ileum[iOrd,c("avg_log2FC","p_val","p_val_adj")]
colnames(df) <- paste0('I-',c('log2FC','p','FDR'))  
df <- cbind(df,DEGs_Colon[iOrd,c("avg_log2FC","p_val","p_val_adj")])
colnames(df)[4:6] <- paste0('C-',c('log2FC','p','FDR'))  
df$gene <- rownames(df)

df$label <- NA
iOrd1 <- which(df$`I-log2FC`> 0.25 & df$`C-log2FC`>0.25) 
iOrd2 <- which(df$`I-log2FC`< -0.25 & df$`C-log2FC` < -0.25)
df$label[c(iOrd1,iOrd2)] <- rownames(df)[c(iOrd1,iOrd2)]

#plot for immune genes:
library(ggrepel)
{
  df$avgFC <- rowMeans(df[,c("I-log2FC","C-log2FC")])
  df<-df[order(df$avgFC),]
  colfunction <- colorRampPalette(c("#3288bd","#e0f3f8","#ffffbf","#fee090","#f46d43","#a50026"))
  df$nodeCol <- colfunction(nrow(df))
  df <- df[df$gene != "CD160",]# remove outliers
  
  p <- ggplot(data=df, aes(x = `I-log2FC`, y = `C-log2FC`,label=label)) +
    geom_point(stroke = 0, alpha = 1,shape = 16,size=3,col = df$nodeCol) + 
    theme_classic() +
    geom_text_repel(colour="black") +
    scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    geom_vline(xintercept=0, linetype="dashed", color = "red")+
    scale_y_continuous(breaks = seq(-6,8,by=0.5))+
    scale_x_continuous(breaks = seq(-4,8,by=0.5))
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
