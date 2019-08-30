# Assembled by Hanna F. Berg for Drug Farm based on and inspired by tutorials from Satijalabs Seurat. 

# A script for loading two samples, merging tem, removing batch effects, clustering and analyzing the differences.

######################################### Load libraries ########################################
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
devtools::load_all("~/dftoolz")

####################################### Load samples ##############################################

ctrl_samp.name="s_SH_s_SHWH_B_cell:Plasma_cell_filt"
stim_samp.name = "u_CHFH_CJL_KUXP_B_cell:Plasma_cell_filt"

ctrl.data<-DGE_load(samp.name = ctrl_samp.name)
stim.data<-DGE_load(samp.name = stim_samp.name)

################# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "successful")
ctrl$stim <- "successful"

# Preform quality control, subset if needed
ctrl[["percent.mt"]] <- PercentageFeatureSet(object = ctrl, pattern = "^mt-")
VlnPlot(object = ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#ctrl <- subset(x = ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 800 & nCount_RNA>500 & nCount_RNA< 1500 & percent.mt < 3.2)
VlnPlot(object = ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Check that ctrl object is scaled and normalized and filtered. If not, do this.

#ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features 
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)


################## Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "unsuccessful", min.cells = 5)
stim$stim <- "unsuccessful"

# Preform quality control, subset if needed
stim[["percent.mt"]] <- PercentageFeatureSet(object = stim, pattern = "^mt-")
VlnPlot(object = stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#stim <- subset(x = stim, subset = nFeature_RNA > 500 & nFeature_RNA <800 & nCount_RNA>500 & nCount_RNA< 1500 & percent.mt < 3 )
VlnPlot(object = stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Check if stim object is scaled and normalized adn filtered. If not, do this.
#stim <- NormalizeData(object = stim, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable features
stim <- FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000)




######################################## Integrate samples #######################################

#Find gene markers to base the integration on
seurat.anchors <- FindIntegrationAnchors(object.list = c(ctrl,stim), dims = 1:30)
seurat.combined <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)

###########
########### Start here for one sample instead of two 
###########
###########

# Add metadata information on which samples are ctrl and stim or successful and unsuccessful.

inf<-ctrl@meta.data[["stim"]]
inf<- append(inf, stim@meta.data[["stim"]])
seurat.combined$stim<-inf

# Or do it like this

seurat.combined$stim <-read.table("path/to/meta-data/file.txt", head =T, row.names = 1)[,1]


############################### Analysis of the integrated sample ##############################

# Re-set the assay to be investigated. We will cluster all cells based on the newly found anchor genes rather than all genes.
DefaultAssay(object = seurat.combined) <- "integrated"

# Find variable genes
seurat.combined <- FindVariableFeatures(object = seurat.combined, selection.method = "vst", nfeatures = 2000)

############# scale data is necessary for every assya in the seurat object, but do it only ONCE per object!
seurat.combined<-ScaleData(seurat.combined)

############## Cluster and visualize

# PCA. npcs can be changed
seurat.combined <- RunPCA(object = seurat.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering, dims and resolution can be changed
seurat.combined <- RunTSNE(object = seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindNeighbors(object = seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)

# visualize clusters. Select which meta data you want to visualize in "group.by" or "split.by".
DimPlot(object = seurat.combined, reduction = "tsne", group.by = "samp.name")
DimPlot(object = seurat.combined, reduction = "tsne", group.by = "stim")
DimPlot(object = seurat.combined, reduction = "tsne", label = TRUE)

DimPlot(object = seurat.combined, reduction = "tsne", split.by = "stim")


################################ Without batch effect removal ###################################

#Test what the samples would look like without the batch effect removal. This creates a new seurat object and does not overwrite all work that has been done.
# also, this is not necessary for the final seurat object, but can be a good comparison to see what the samples would look like with batch effects.

assay.RNA<-GetAssay(object=seurat.combined, assay="RNA")
DGE.wbatch<-as.matrix(GetAssayData(object=assay.RNA, slot = "data"))
seurat.combined.wbatch<-CreateSeuratObject(counts=DGE.wbatch, project="with batch effects")
seurat.combined.wbatch$stim<-seurat.combined@meta.data[["stim"]]

seurat.combined.wbatch <- FindVariableFeatures(object = seurat.combined.wbatch, selection.method = "vst")
seurat.combined.wbatch<-ScaleData(seurat.combined.wbatch)

#make sure these parameters are the same as in the object with the reduced batch effects, otherwise they will not be comparable.
seurat.combined.wbatch <- RunPCA(object = seurat.combined.wbatch, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
seurat.combined.wbatch <- RunTSNE(object = seurat.combined.wbatch, reduction = "pca", dims = 1:30)
seurat.combined.wbatch <- FindNeighbors(object = seurat.combined.wbatch, reduction = "pca", dims = 1:30)
seurat.combined.wbatch <- FindClusters(seurat.combined.wbatch, resolution = 1)

# visualize clusters
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", group.by = "samp.name")
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", group.by = "stim")
DimPlot(object = seurat.combined.wbatch, reduction = "tsne", label = TRUE)


DimPlot(object = seurat.combined.wbatch, reduction = "tsne", split.by = "stim")

# You are now done with the comparison of the samples with batch effects. The rest of the script analyses the integrated, batch effect removed seurat object.

##################################### Differential expression ctrl vs stim #################################

##############If you already computed and saved the DE genes, load it now and skip to the Visualization part.
seurat.combined.markers<-read.table("/location/for/Gene markers samples.txt", head=T, sep="\t")

############### If this is the first time computing the DE genes, continue below:

# Re-set the assay to be investigated. We will look at differential expression in all genes rather than the anchor genes.
DefaultAssay(object = seurat.combined) <- "RNA"

# Choose what differential genes you want to compute. Between clusters, between successful vs unsuccessful, between samples and so on.
Idents(object = seurat.combined)<- as.factor(seurat.combined@meta.data[["samp.name"]]) 

# DO ONLY ONCE: scale ALL RNA data to make sure to include all genes when searching for DE genes. 
seurat.combined <- ScaleData(seurat.combined,features = rownames(GetAssay(seurat.combined, slot = "RNA")))


# Find all markers calculates the differential expression. This takes a long time for big data sets, make sure you save the DE genes after this.
seurat.combined.markers <- FindAllMarkers(seurat.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Save DE genes
write.table(seurat.combined.markers, file="Gene markers samples.txt", sep ="\t")


############################## Visualize the differentially expressed genes ####################################

#Select how many genes that will be plotted. 
top10 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# or:
top50 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

############## Heatmap

#the heatmap can not plot more than about 18000 cells in one go. Create subsets if the seurat.combined object contains more cells than 18000

seurat_small<-subset(x=seurat.combined, downsample=8000 )

pdf("Heatmap.pdf")
DoHeatmap(seurat_small, features=as.character(top50$gene), size = 2) + theme(axis.text.y = element_text(size = 1))
dev.off()

# if seurat.combined dontains less than 18000 cells, plot all cells.
DoHeatmap(seurat.combined, features = as.character(top50$gene), slot = "scale.data", size = 3) + theme(axis.text.y = element_text(size = 5))+ggtitle(paste("Heatmap", ctrl_samp.name, "vs", stim_samp.name))

############### Create new meta.data
# To compare for example all successful and all unsuccessful within one cluster, this data has to be created. 
# The function below assigns a name such as "cluster 1, successful" "and cluster 1 unsuccessful" for each cell. Then re-run findAllMarkers and heatmap to visualize.

samp.clus.ident=c()

for (i in 1:length(seurat.combined@active.ident)){
  nw<-paste0(seurat.combined@meta.data[["seurat_clusters"]][i]," ",seurat.combined@meta.data[["stim"]][i])
  samp.clus.ident<-append(samp.clus.ident, nw, after = i)
}

# set the new meta data to the identity to be evaluated.
Idents(object = seurat.combined)<- as.factor(samp.clus.ident)

############### Heatmap of individual clusters
#Plot a heatmap of only cluster x to compare successful and unsuccessful.

# Create an empty vector only once.
top_DE_ctrl_stim<-data.frame()

# subset the clusters to compare
seurat_subset<-subset(x = seurat.combined, idents = c("1 successful", "1 unsuccessful"))

# Find differentially expressed genes.
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

# Plot heatmap
pdf("Heatmap 1", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

# extract the most differentiated genes between successful and unsuccessful in this cluster.
top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

# save those genes in a data frame.
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


############################################### Dot plot #############################################
# Dot plot is a quick and useful tool to compare some specific genes between samples, clusters or whatever is in the active identity.

# Plot top differential genes
DotPlot(object = seurat.combined, features = rev(x = unique(top10$gene)), cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))

# Plot the markers with the biggest difference between successful and unsuccessful within one cluster. 
DotPlot(object = seurat.combined, features = rev(x = unique(top_DE_ctrl_stim$gene)), cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))

#plot selected features and title
features = c("GENE1", "GENE2", "GENE3")
DotPlot(object = seurat.combined, features = features, cols = c("blue", "red"),  dot.scale = 8, split.by = "stim") + theme(axis.text.x = element_text(size = 7, angle = 90))+ggtitle(paste("DotPlot", ctrl_samp.name, "vs", stim_samp.name))

# Add a vertical line
DotPlot(seurat.combined,features = features,cex.use=4)+ theme(axis.text.x = element_text(size = 7, angle = 90))+ geom_vline(xintercept=8.5, color="red")

############################################## Save seurat #########################################
samp.name= "The name of my sample"
filen<-paste0(samp.name,".Seurat.rds")
saveRDS(seurat.combined, file=filen)


