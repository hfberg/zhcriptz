# a script to filter and normalize DGE:s before merging. 
#The DGEs produced here should alway be called it's name_filt. They have to have sample names in their barcode and be ready for merging.

devtools::load_all("~/dftoolz")

library(Seurat)
library(dplyr)
samp.name = "u_CJL_pre1_B_cell:Plasma_cell"

#if DGE doesn't have samp names, set  name.barc=T in DGE_load
raw.data<- DGE_load("/home/drugfarm/proj.zhqugen/BM/process (copy)/data/u_CJL_pre1_plasma_dge.txt",path = T)


seurat <- CreateSeuratObject(counts = raw.data, project = samp.name)

seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat <- subset(x = seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 4800 & nCount_RNA>1000 & nCount_RNA< 18000 &percent.mt < 10)

#seurat <- subset(x = seurat, subset = nFeature_RNA > 1000  & nCount_RNA>1000 &percent.mt < 10)


VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
DGE<-as.matrix(GetAssayData(object = seurat))
samp.name.filt<-paste0(samp.name,"_filt")
DGE_save(DGE, samp.name=samp.name.filt)
print(samp.name.filt)