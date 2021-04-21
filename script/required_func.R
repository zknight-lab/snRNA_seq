###initialize object###
library(Seurat)
library(glmGamPoi)

Rm_dblet <- function(counts, path){
  seurat <- init_Seurat(counts = counts)
  scrblt <- read.table(paste0(path, 'doublet_prediction.txt'), head = TRUE, sep = ',')
  Doublet <- subset(scrblt, scrblt$prediction == 'True')
  counts <- seurat@assays$RNA@counts
  rname <- rownames(counts)
  cname <- colnames(counts)
  ##nuclear seq
  counts <- counts[-which(grepl(pattern = '^mt-', rname) == TRUE), -which(cname %in% Doublet$barcodes)]
  #counts <- counts[, -which(cname %in% Doublet$barcodes)])
  return(counts)
}

init_Seurat <- function(counts, min.cells = 3, min.features = 3){
  seurat <- CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features)
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  return(seurat)
}

Create_Cluster <- function(seurat){
  seurat <- NormalizeData(seurat, normalization.method = normalization.method , scale.factor = scale.factor)
  seurat <- FindVariableFeatures(seurat, selection.method = selection.method, nfeatures = nfeatures)
  seurat <- ScaleData(seurat,features=VariableFeatures(seurat)) 
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  seurat <- FindNeighbors(seurat, dims = 1:dim)
  seurat <- FindClusters(seurat, resolution = 0.5)
  seurat <- RunTSNE(seurat, reduction="pca",dims = 1:dim,seed.use = seed.use)
  seurat <- RunUMAP(seurat, dims = 1:dim)
  return(seurat)
}

Create_SCT_Cluster <- function(seurat){
  seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
  seurat <- SCTransform(seurat,  vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = FALSE)
  seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
  seurat <- FindNeighbors(seurat, dims = 1:dim)
  seurat <- FindClusters(seurat, resolution = 0.5)
  seurat <- RunTSNE(seurat, reduction="pca",dims = 1:dim,seed.use = seed.use)
  seurat <- RunUMAP(seurat, dims = 1:dim)
  return(seurat)
}

Create_maker <- function(seurat, top_n = 20, only.pos = TRUE, min.pct = 0.25){
  all.markers <- FindAllMarkers(seurat, only.pos = only.pos, min.pct = min.pct, logfc.threshold = 0.25)
  return(all.markers)
}

###################
###Plot Function###
###################
library(ggplot2)
library(ggsci)
library(dplyr)

plot_QC <- function(seurat){
  plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2  
  VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
}

plot_PCA <- function(seurat){
  ElbowPlot(seurat, ndims = 20,reduction = "pca")
  
}

plot_TSNE <- function(tsne, labels, color, size = 0.75, alpha = 1){
  cbPalette <- c("#305BA9", "#EF7F37", "#BEDA7F", "#BB002199", "#8365A9", "#FDA8A7", "#E69F00", "#FF2121", "#63D8B7", "#00C4FF", "#63187999", "#9999FF", "#B84664","#66CCFF","#FCB045","#EE000099","#EC692A","#1ABC9C", "#2ECC71", "#3498DB", "#9B59B6", "#34495E", "#16A085", "#27AE60", "#2980B9", "#8E44AD", "#2C3E50")
  cbPalette_1 <- c("#1ABC9C", "#2ECC71", "#3498DB", "#9B59B6", "#34495E", "#16A085", "#27AE60", "#2980B9", "#8E44AD", "#2C3E50", "#F1C40F", "#E67E22", "#E74C3C", "#ECF0F1", "#95A5A6", "#F39C12", "#D35400", "#C0392B", "#BDC3C7", "#7F8C8D")
  ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = color)) + 
    geom_point(size = size, alpha = alpha) + 
    scale_color_manual(labels = labels, values = cbPalette) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
    ) +
    guides(color=guide_legend(title = "subgroup", override.aes=list(size=5)))
}


plot_PDF <- function(seurat, filename){
  tsne = seurat@reductions$tsne@cell.embeddings %>%
    as.data.frame() %>% cbind(tx = seurat@meta.data$seurat_clusters,ident = seurat@meta.data$orig.ident)
  DimPlot(seurat, reduction = "umap")+ggtitle("UWOT")+theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(filename, '_umap.pdf'), device = pdf, path = output, width = 8.3, height = 8.3)
  
  plot_TSNE(tsne, sort(unique(tsne$tx)),tsne$tx)
  ggsave(filename = paste0(filename, '_tsne_clusters.pdf'), device = pdf, path = output, width = 8.3, height = 8.3)
  
  VizDimLoadings(seurat, dims = 1:10, reduction = "pca",nfeatures = 10)
  ggsave(filename = paste0(filename, '_pca.pdf'), device = pdf, path = output, width = 13.3, height = 8.3)
  
  DimPlot(seurat, reduction = "pca")
  ggsave(filename = paste0(filename, '_pcadimplot.pdf'), device = pdf, path = output, width = 6.3, height = 6.3)
  
  (DimHeatmap(seurat, dims = 1:6, cells = 500, balanced = TRUE)) %>%
    ggsave(filename = paste0(filename, '_Heatmap.pdf'), device = pdf, path = output, width = 8.3, height = 8.3)
}

plot_Feature <- function(seurat, all.marker, top_n = 20, path, filename){
  dge<-all.markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC) 
  write.csv(dge, path = paste0(path, filename, '_DGE.csv'))
  
  top5 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  size = length(unique(seurat$seurat_cluster))
  FeaturePlot(seurat, features = top5$gene, reduction = "tsne", ncol = 5)
  ggsave(filename = paste0(filename, '_FeaturePlot.pdf'), device = pdf, path = output, width = 4*5, height = 4*size, limitsize = FALSE)
  
  top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  DoHeatmap(seurat, features = top10$gene, label = FALSE) 
  ggsave(filename = paste0(filename, '_Heatmap_top10$gene.pdf'), device = pdf, path = output, width = 1.6*size, height = 0.9*size, limitsize = FALSE)
}