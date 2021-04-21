rm(list = ls())

source('d:/snRNA_seq/script/required_func.R')
source('d:/snRNA_seq/script/config.R')

counts <- Read10X(data.dir = paste0(datapath, filename))
output <- paste0(output, filename, '/')
counts <- Rm_dblet(counts, output)

brain_si <- init_Seurat(counts)

plot_QC(brain_si)

brain_si <- Create_SCT_Cluster(brain_si)

plot_PDF(brain_si, paste0(filename, '_si'))

all.marker <- Create_maker(brain_si)
#cluster1.markers <- FindMarkers(brain_si, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc")
##----
plot_Feature(brain_si, all.marker)
save(list=ls(),file = paste0(RData,filename, "_si.RData"))


marker <- c( 'Map2', 'Mbp', 'Slc17a6', 'Slc32a1', 'Chat', 'Slc18a3','Ddc', 'Slc6a3', 'Slc18a2')

x <- brain_si$seurat_clusters

x <- factor(x, levels = sort(unique(x), na.last = TRUE),labels = c('Neuron_0', 'Neuron_1', 'Neuron_2', 'Oligodendrocytes', 'Astrocytes', 'Neuron_3', 'Neuron_4', 'Neuron_5', 'Microglia', 'Neuron_6', 'Endothelial Cells_1', 'Type IC spiral ganglion neuron', 'Cholinergic neuron', 'Endothelial Cells_2'), exclude = NA, ordered = is.ordered(x))

brain_si$seurat_clusters <- x

top <- brain_si.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
VlnPlot(brain_si, features = c(top$gene),stack=T,same.y.lims = T)+theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + NoLegend()+ylab(label = 'Cluster')
ggsave(filename = 'si_without_mt_VlnPlot.pdf', device = pdf, path = output, width = 24.9, height = 8.3, limitsize = FALSE)

VP<-VlnPlot(brain_si, features =marker,stack=T,same.y.lims = T,fill.by = 'ident',pt.size = 0.5)+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),element_line(linetype = 'blank'),panel.border = element_blank())+NoLegend()+scale_fill_brewer(palette = 'Set1')
VP$plot_env$geom[[1]]$geom$default_aes$linetype<-c('blank')
VP<-VP+theme(strip.text.x = element_text(angle = 75))
VP

#fretch Seurat data
exprs <- data.frame(FetchData(object = brain_si, vars = gene))
exprs$Barcode <- rownames(exprs)
clust<-data.frame()
#barcode与聚类信息提取
ident<-data.frame(Barcod=names(newname.exprs$Barcode),orig.ident=newname.brain_si$seurat_clusters)
#通过merge函数，将表达量与聚类号对应起来
c<-merge(exprs,ident,by='Barcod')
#对其进行排序
c$orig.ident<-factor(c$orig.ident,levels=c(sort(unique(immune.combined@ident))))


ggplot(data,aes(x=factor(orig.ident),y=Cdh17,fill=orig.ident))+geom_violin(alpha=0.8,width=1)
ggsave('test.pdf')

noise <- rnorm(n = length(x = data[,c('Cdh17')])) / 100000
data[,c('Cdh17')] <- data[, c('Cdh17')] + noise
ggplot(data = data,mapping = aes(x = factor(x = orig.ident),y = Cdh17)) +geom_violin(scale = "width",adjust =1,trim = TRUE,mapping = aes(fill = factor(x = orig.ident))) 
ggsave('test.pdf')

