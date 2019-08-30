# A script to create different types of heatmaps on differentially expressed genes in clusters.

setwd("/home/drugfarm/proj.zhqugen/BM/process")
devtools::load_all("~/dftoolz")
library(SingleR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

raw<- DGE_load("s_seurat_SH_SHWH+u_CHFH_CJL_KUXP_B_cell:Plasma_cell", name.barc = F)
annot<-annot_load(samp.name="s_seurat_SH_SHWH+u_CHFH_CJL_KUXP_B_cell:Plasma_cell_filtered", cluster = F,samp = T)

seurat<-CreateSeuratObject(raw, meta.data=annot)
seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat <- subset(x = seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & nCount_RNA>1000 & nCount_RNA< 12500 & percent.mt < 10 )
VlnPlot(object = seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat<- NormalizeData(seurat)

seurat <- ScaleData(seurat, features = rownames(GetAssayData(seurat, slot = "data")))

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)



# Find top genes, extract to dataset
Idents(object = seurat) <- annot
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top30 <- seurat.combined.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
top50 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

# subset seurat to plot in heatmap.
seurat_small<-subset(x=seurat.combined, downsample=1000 )
pdf("PDF1_view_before_deleate.pdf")
DoHeatmap(seurat_small, features=as.character(top30$gene), size = 2) + theme(axis.text.y = element_text(size = 2))
dev.off()

#compare successful to unsuccessful

annot_succ<-read.csv("/home/drugfarm/proj.zhqugen/PBMC/process/data/all_T_cell:gamma-delta.SUCC.UNSUCC.txt", row.names = 1, header = T)

#repeat above steps with idents set to annot_succ instead.
Idents(object = seurat) <- annot_succ

####################### Heatmap of all clusters #####################################################

top_DE_ctrl_stim<-data.frame()

seurat_subset<-subset(x = seurat.combined, idents = c("1 successful", "1 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 1", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))

seurat_subset<-subset(x = seurat.combined, idents = c("2 successful", "2 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 2", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("3 successful", "3 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 3", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("4 successful", "4 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 4", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("5 successful", "5 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 5", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("6 successful", "6 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 6", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("7 successful", "7 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)



pdf("Heatmap 7", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("8 successful", "8 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 8", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("9 successful", "9 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 9", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("10 successful", "10 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 10", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("11 successful", "11 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 11", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("12 successful", "12 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 12", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("13 successful", "13 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 13", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("14 successful", "14 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 14", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("15 successful", "15 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 15", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("16 successful", "16 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 16", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("17 successful", "17 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 17", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("18 successful", "18 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 18", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("19 successful", "19 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 19", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("20 successful", "20 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 20", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("21 successful", "21 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 21", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("22 successful", "22 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 22", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("23 successful", "23 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 23", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))


seurat_subset<-subset(x = seurat.combined, idents = c("24 successful", "24 unsuccessful"))
markers_subset <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25)

pdf("Heatmap 24", height = 8.50, width = 11)
DoHeatmap(seurat_subset, features = markers_subset$gene, slot = "scale.data", size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 3))
dev.off()

top1_sub <- markers_subset %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
top_DE_ctrl_stim<- rbind(top_DE_ctrl_stim, as.data.frame(top1_sub))

####################################################################################################
