library(SingleR)
library(Seurat)
library(ggplot2)
library(scater)
library(dplyr)
library(pheatmap)
library(ggsci)

rm(list = ls())
source('d:/snRNA_seq/script/config.R')
source('d:/snRNA_seq/script/required_func.R')
load(paste0(RData,'brain_ref.RData'))
load(paste0(RData,'brain_si_tdTomato.RData'))

output = paste0(output, 'tdTomato')
for (i in 0:22) {
  test = brain_si@assays$RNA@counts[,which(brain_si$seurat_clusters == i)]
  pred_test <- SingleR(test = test, ref = ref, aggr.ref = FALSE, labels = ref$label.main)
  print(paste('This is ', i, ' cluster.'))
  print('-------------------------------')
  (plotScoreHeatmap(pred_test, annotation_col = as.data.frame(brain_si$seurat_clusters))) %>%
      ggsave(filename = paste0('score_ref-Mou_cluster', i, '.pdf'), device = pdf, path = output, width = 6.3, height = 9.3)
  #ggsave(filename = paste0('score_cluster_', i, '.pdf'), device = pdf, path = output, )
  print(table(pred_test$labels))
  print(summary(pred_test$labels))
  # png(paste0(workpath, output, '\\', 'cluster_', i, '_score_distribution.png'), res = 300, width = 3000, height = 2500)
  plotScoreDistribution(pred_test, ncol = 3, show.nmads = 3)
  ggsave(filename = paste0('distribution_ref-Mou_cluster', i, '.pdf'), device = pdf, path = output, width = 6.3, height = 21.3)
}

test = brain_si@assays$RNA@counts
pred_test <- SingleR(test = test, ref = ref_Mou, aggr.ref = FALSE, labels = ref_Mou$label.main)

(plotScoreHeatmap(pred_test, annotation_col = as.data.frame(brain_si$seurat_clusters)) ) %>%
  ggsave(filename = paste0('cluster_label', '.pdf'), device = pdf, path = output, width = 9, height = 9)

plotScoreDistribution(pred_test, ncol = 3)

all.markers <- metadata(pred_test)$de.genes
test@factor <- pred_test$labels

plotHeatmap(test, order_columns_by = 'labels', features = unique(unlist(all.markers)))

x <- brain_si$seurat_clusters

x <- factor(x, levels = sort(unique(x), na.last = TRUE),labels = c('Astrocytes_1', 'Neuron_1', 'Neuron_2', 'Oligodendrocytes', 'Astrocytes_2', 'Neuron_3', 'Neuron_4', 'Neuron_5', 'Microglia', 'Neuron_6', 'Endothelial Cells_1', 'Type IC spiral ganglion neuron', 'Cholinergic neuron', 'Endothelial Cells_2'), exclude = NA, ordered = is.ordered(x))

brain_si$seurat_clusters <- x

labels = c('Neuron_0', 'Neuron_1', 'spiral ganglion neuron', 'Oligodendrocytes', 'Astrocytes', 'Neuron_3', 'Neuron_4', 'Neuron_5', 'Microglia', 'Neuron_6', 'Endothelial Cells_1', 'Neuron_7', 'Neuron_8', 'Endothelial Cells_2')
tsne = brain_si@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = pred_test$labels,cluster = brain_si@meta.data$seurat_clusters)


source('d:/snRNA_seq/plot_func.R')

non_neuron <- c(0,2,13,14,19)
tsne_non <- tsne[which(tsne[,4] %in% non_neuron),]

tsne_neu <- tsne[which(! tsne[,4] %in% non_neuron),]
plot_TSNE(tsne, sort(unique(tsne$tx)), tsne$tx) +
plot_TSNE(tsne_neu, sort(unique(tsne_neu$tx)), tsne_neu$tx)


neuron@meta.data$seurat_clusters <- as.factor(pred_test$labels)
counts <- neuron@assays$RNA@counts[, which(neuron@meta.data$seurat_clusters == 'Neurons')]

counts <- neuron@assays$RNA@counts[,which(! neuron@meta.data$seurat_clusters %in% non_neuron)]

neuron <- init_Seurat(counts, min.cell = 3, min.features = 3)
plot_QC(neuron)
neuron <- Create_SCT_Cluster(neuron)
filename = 'tdTomato'
plot_PDF(neuron, paste0(filename, '_neuron'))

marker = c('Slc17a6', 'Slc32a1', 'Map2', 'Pvalb', 'Ddc', 'Slc18a2', 'Gad1', 'Gad2', 'Chat', 'tdTomato')
FeaturePlot(neuron, features = marker, reduction = 'tsne', ncol = 2)

Chat = brain_si

filename = 'chat_Cluster'
all.markers <- FindAllMarkers(CHAT_, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dge<-all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
##----
top <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

FeaturePlot(CHAT_, features = top$gene,reduction = "tsne", ncol = 5)
ggsave(filename = paste0(filename, '_SCT_FeaturePlot.pdf'), device = pdf, path = output, width = 16, height = 9, limitsize = FALSE)

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(CHAT_, features = top10$gene, label = FALSE) 
ggsave(filename = paste0(filename, '_SCT_Heatmap_top10$gene.pdf'), device = pdf, path = output, width = 9, height = 9)

save(list=c('all.markers', 'dge', 'CHAT_'),file = paste0(RData,filename,'_tdTomato.RData'))

##calculate averageExp

AverageExp<-AverageExpression(CHAT_,features=unique(dge$gene))
head(AverageExp$RNA)
library(psych)
library(pheatmap)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)