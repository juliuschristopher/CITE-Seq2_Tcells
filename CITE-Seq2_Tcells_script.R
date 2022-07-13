## ##CITE-Seq2 T cells - CD4+ and CD8+ T cells####
#Julius Christopher Baeck

####Set Up####
setwd("/Volumes/GoogleDrive/Shared drives/Okkengroup/Experiments/Julius/Experiments/CITE-Sequencing/CITE-Seq (2)/Overall_analysis/CITE-Seq2_Tcells/CITE-Seq2_Tcells")

##Packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(RColorBrewer)
library(writexl)
library(ggridges)
library(clustree)
library(scRepertoire)
library(future)
library(alakazam)
library(immunarch)
library(airr)
library(biomaRt)
library(SeuratDisk)
library(SeuratData)
library(stringr)
library(viridis)
library(escape)
library(dittoSeq)
library(SingleCellExperiment)
library(ggsci)
library(pals)
library(harmony)
library(gridExtra)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(org.Mm.eg.db)
library(ggpubr)
library(data.table)
library(Polychrome)

##Functions##
tfidf = function(data,target,universe){
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,target,drop=FALSE]>0)
  nTot = Matrix::rowSums(data[,universe,drop=FALSE]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,target,drop=FALSE]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[,universe[!(universe%in%target)],drop=FALSE]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}

##Colours##
col_con <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_con1 <- viridis(50)
col_con2 <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
col_con2 <-as.character(col_con2)

swatch(col_con2)

##Load precessed seurat object - CITESeq2_plain.h5seurat##
CITESeq2_plain.h5seurat <- LoadH5Seurat("CITESeq2_plain.h5seurat") #30843 cells
head(CITESeq2_plain.h5seurat[[]])

Total_cells <- CITESeq2_plain.h5seurat
head(Total_cells[[]])
meta_data <- Total_cells[[]]

summary <- meta_data %>%
  group_by(orig.ident) %>%
  summarize(n())

##Add mice and genotype information##
Idents(Total_cells) <- Total_cells$orig.ident
Total_cells <- RenameIdents(Total_cells, `a` = "WT_1", `b` ="WT_2", `c` ="BCL6_1", `d` = "BCL6_2", `f` = "E1020K_1",`g` ="E1020K_BCL6_1", `h` = "E1020K_BCL6_2")
Total_cells[["Mouse"]] <- Idents(Total_cells)
Idents(Total_cells) <- Total_cells$Mouse

Idents(Total_cells) <- Total_cells$orig.ident
Total_cells <- RenameIdents(Total_cells, `a` = "WT", `b` ="WT", `c` ="BCL6", `d` = "BCL6", `f` = "E1020K",`g` ="E1020K_BCL6", `h` = "E1020K_BCL6")
Total_cells[["Genotype"]] <- Idents(Total_cells)

Idents(Total_cells) <- Total_cells$Mouse

####Filter out T cells###
test <- Total_cells #30843

#Based on no clonotype, ADT CD4 >2, ADT CD8 > 2 and RNA Cd19 < 1 expression
DefaultAssay(test) <-  "ADT"
T_cells <- subset(test, subset = cloneType == "NA") #15651
T_cells <- subset(T_cells, subset = Cd4 > 2 | Cd8a > 2) #14140

DefaultAssay(T_cells) <-  "RNA"
T_cells <- subset(T_cells, subset = Cd19 < 1) #13377


####Dimensionality reduction and clustering####
##Normalise subset##
DefaultAssay(T_cells) <- "RNA" #For log normalisation

#RNA normalisation#
T_cells <- NormalizeData(T_cells, verbose = TRUE)
T_cells <- FindVariableFeatures(T_cells, nfeatures = 3000)
T_cells <- ScaleData(T_cells)

#Visualisation#
top20 <-  head(VariableFeatures(T_cells), 20)
plot1.1 <-  VariableFeaturePlot(T_cells)
top20_plot <-  LabelPoints(plot = plot1.1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
top20_plot

#RNA PCA#
T_cells <- RunPCA(T_cells, verbose = FALSE, features = VariableFeatures(object = T_cells))
pca_variance <- T_cells@reductions$pca@stdev^2
plot(pca_variance/sum(pca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #24

#RNA clustering#
DefaultAssay(T_cells) <- "RNA" #For log normalisation

T_cells <- FindNeighbors(T_cells, dims = 1:24)
T_cells <- FindClusters(T_cells, resolution = 0.3, verbose = FALSE) #0.3 for the resolution
clustree(T_cells, prefix = "RNA_snn_res.") + theme(legend.position="bottom")
T_cells <-RunUMAP(T_cells, dims = 1:24, assay = 'RNA', reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
T_cells_p1 <- DimPlot(T_cells, label = TRUE, reduction = "rna.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("RNA Clustering") + theme_bw() + NoLegend()
T_cells_p1 <- T_cells_p1 + theme(plot.title = element_text(color="black", size=25, face="bold"))

DefaultAssay(T_cells) <- "RNA"
FeaturePlot(T_cells, features = "Gzma", reduction = "rna.umap", shape.by = "Genotype", pt.size = 1)

?FeaturePlot


#ADT#
DefaultAssay(All_cells) <- "ADT"

#ADT normalisation#
VariableFeatures(All_cells) <- rownames(All_cells[["ADT"]])
All_cells <- NormalizeData(All_cells, normalization.method = "CLR", margin = 2)
All_cells <- ScaleData(All_cells)

#ADT PCA#
All_cells <- RunPCA(All_cells, reduction.name = 'apca', approx = FALSE)
apca_variance <- All_cells@reductions$apca@stdev^2
plot(apca_variance/sum(apca_variance), 
     ylab="Proportion of variance explained", 
     xlab="Principal component")
abline(h = 0.01) #24

#ADT clustering
All_cells <- FindNeighbors(All_cells, dims = 1:24, reduction = "apca")
All_cells <- FindClusters(All_cells, resolution = 1.2, verbose = FALSE) #1.2 for the resolution
clustree(All_cells, prefix = "ADT_snn_res.") + theme(legend.position="bottom")
All_cells <- RunUMAP(All_cells, reduction = 'apca', dims = 1:24, assay = 'ADT', reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
All_cells_p2 <- DimPlot(All_cells, label = TRUE, reduction = "adt.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("ADT Clustering") + theme_bw() + NoLegend()
All_cells_p2 <- All_cells_p2 + theme(plot.title = element_text(color="black", size=25, face="bold"))

#WNN#
DefaultAssay(All_cells) <- "RNA" #For log normalisation

#Combine into wnn plot#
All_cells <- FindMultiModalNeighbors(
  All_cells, reduction.list = list("pca", "apca"), 
  dims.list = list(1:25, 1:24), modality.weight.name = "RNA.weight")

#WNN clustering#
All_cells <- FindClusters(All_cells, graph.name = "wsnn", algorithm = 3, resolution = 0.7, verbose = TRUE) #0.7 for the resolution
clustree(All_cells, prefix = "wsnn_res.") + theme(legend.position="bottom")
All_cells <- RunUMAP(All_cells, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
All_cells_p3 <- DimPlot(All_cells, label = TRUE, reduction = "wnn.umap", pt.size = 1.3, label.size = 6, label.box = TRUE, cols = col_con2) +  ggtitle("Seurat Clusters") + theme_bw() + NoLegend()
All_cells_p3 <- All_cells_p3 + theme(plot.title = element_text(color="black", size=25, face = "bold")) + xlab("UMAP1") + ylab("UMAP2")

All_cells_p1 + All_cells_p2 + All_cells_p3

####No batch effect due to all samples coming from the same batch####
head(All_cells[[]])

#Save All_cells Seurat object#
SaveH5Seurat(All_cells, filename = "CITE-Seq2_all_cells", overwrite = TRUE)
