options(connectionsObserver = NULL)
library(RSQLite)
library(clusterProfiler)
library(topGO)
library(org.Mm.eg.db)
library(Rgraphviz)
library(pathview)
library(enrichplot)
library(DOSE)
library(ggnewscale)
library(ggsci)

id_to_symbol <- function(result){
  for (i in 1:length(result)){
    id <- unlist(strsplit(result[i],"/"))
    symbol <-mapIds(x = org.Mm.eg.db,
                    keys = id,
                    keytype = "ENTREZID",
                    column = "SYMBOL")
    symbol <- paste(symbol,collapse = "/")
    result[i] <- symbol
  }
  return(result)
}


load('D:/snRNA_seq/RData/chat_Cluster_tdTomato.RData')

##symbol to entrez_id
Sym2Ent <- function(dge){
  symbol <- read.table('d:/snRNA_seq/symbol.txt', sep = ',', header = TRUE)[,-1]
  DGE.ensembl_id <- symbol[which(symbol$sym %in% dge$gene),1]
  DGE.entrez_id <- mapIds(x = org.Mm.eg.db,
                          keys = as.character(DGE.ensembl_id),
                          keytype = "ENSEMBL",
                          column = "ENTREZID")
  DGE <- dge %>% cbind(ENSEMBL = DGE.ensembl_id, ENTREZID = DGE.entrez_id)
  DGE <- DGE[which(! is.na(DGE$ENTREZID)),]
  #DEG.entrez_id = na.omit(DEG.entrez_id)
  return(DGE)
}

DGE <- Sym2Ent(dge)

gene_erichment_results = list()
universe <- CHAT_@assays$RNA@counts@Dimnames[1]
for (cl in as.character(unique(levels(DGE$cluster)))){
  entrez_id <- subset(DGE, cluster == cl)$ENTREZID
  logFC <- subset(DGE, cluster == cl)$avg_log2FC
  # gene_erichment_results[[cl]] = list()
  for (ont in c('BP', 'MF', 'CC')) {
    erich.go <- enrichGO(gene = entrez_id,
                            OrgDb = org.Mm.eg.db,
                            keyType = "ENTREZID",
                            ont = ont,
                            pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
    gene_erichment_results[[cl]][[ont]] = erich.go
  
    # dotplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))
    # ggsave(paste0('d:/snRNA_seq/output/tdTomato/',ont, '_chat_dotplot_', cl, '.pdf'), device = pdf,width = 9, height = 9)
    # barplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))
    # ggsave(paste0('d:/snRNA_seq/output/tdTomato/',ont, '_chat_barplot_', cl, '.pdf'), device = pdf,width = 9, height = 9)
    # 
    # cnetplot(gene_erichment_results[[cl]][[ont]], categorySize="pvalue", showCategory = 10, foldChange=logFC, colorEdge = TRUE)
    # ggsave(paste0('d:/snRNA_seq/output/tdTomato/',ont, '_chat_cnet_', cl, '.pdf'), device = pdf,width = 9, height = 9)

  }
}
answer <- gene_erichment_results

for (cl in as.character(unique(levels(DGE$cluster)))){
  for (ont in c('BP', 'CC')) {
    dotplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))
    ggsave(paste0('d:/snRNA_seq/output/tdTomato/',cl, '_chat_dotplot_', ont, '.pdf'), device = pdf,width = 18, height = 12)
    barplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))
    ggsave(paste0('d:/snRNA_seq/output/tdTomato/',cl, '_chat_barplot_', ont, '.pdf'), device = pdf,width = 18, height = 12)

    cnetplot(gene_erichment_results[[cl]][[ont]], categorySize="pvalue", showCategory = 10, foldChange=logFC, colorEdge = TRUE) + theme(text = element_text(size = 12, face = "bold"))
    ggsave(paste0('d:/snRNA_seq/output/tdTomato/',cl, '_chat_cnet_', ont, '.pdf'), device = pdf,width = 18, height = 12)
    cnetplot(gene_erichment_results[[cl]][[ont]], circular = TRUE, categorySize="pvalue", showCategory = 10, foldChange=logFC, colorEdge = TRUE) +theme(text = element_text(size = 12, face = "bold"))
    ggsave(paste0('d:/snRNA_seq/output/tdTomato/',cl, '_chat_circle_', ont, '.pdf'), device = pdf,width = 18, height = 12)
    
  }
}
 
library(cowplot)

for (cl in as.character(unique(levels(DGE$cluster)))){
  for (ont in c('BP', 'CC')) {
    p1 <- dotplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))
    p2 <- barplot(gene_erichment_results[[cl]][[ont]]) + theme(axis.text.y = element_text(size = 12, face = "bold"))

    p3 <- cnetplot(gene_erichment_results[[cl]][[ont]], categorySize="pvalue", showCategory = 10, foldChange=logFC, colorEdge = TRUE) + theme(text = element_text(size = 12, face = "bold"))
    p4 <- cnetplot(gene_erichment_results[[cl]][[ont]], circular = TRUE, categorySize="pvalue", showCategory = 10, foldChange=logFC, colorEdge = TRUE) +theme(text = element_text(size = 12, face = "bold"))
    gg <- ggdraw() +     
      draw_plot(p1, x=0, y=0.5, width=0.5, height=0.5) +  
      draw_plot(p2, 0, 0, 0.5, 0.5) +   
      draw_plot(p3, 0.5, 0, 0.5, 0.5) + 
      draw_plot(p4,0.5, 0.5, 0.5, 0.5)
    gg
    ggsave(paste0('d:/snRNA_seq/output/tdTomato/',cl, '_chat_total_', ont, '.pdf'), device = pdf,width = 28, height = 21)
    
  }

}

## 画图
barplot(erich.go.CC) + theme_gray() + scale_color_gradient(low = 'green', high = 'red')

#plotGOgraph(erich.go.BP)
pdf(file="e:/Cleandata/enrich.go.bp.tree.pdf",width = 10,height = 15)
plotGOgraph(gene_erichment_results[[1]][['CC']])
dev.off()

goplot()
erich.kegg <- enrichKEGG(gene = DEG.entrez_id,
                       organism = "mmu",
                       keyType = "kegg",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

re <- id_to_symbol(erich.kegg@result$geneID) 
erich.kegg@result$geneID <- re

re <- id_to_symbol(erich.go.BP@result$geneID) 
erich.go.BP@result$geneID <- re

barplot(erich.kegg)
##横轴为该pathway的差异基因个数，
##纵轴为富集到的pathway的描述信息， 

cnetplot(erich.kegg,showCategory = 10)
##showCategory指定展示的pathway的个数，默认展示显著富集的top10个，即p.adjust最小的10个。
##图中点的颜色对应p.adjust的值，从小到大，对应蓝色到红色，大小对应该GO terms下的差异基因个数，个数越多，点越大。

cnetplot(erich.kegg,showCategory = 10,circular=T,  ###画为圈图
         colorEdge=T)      ##线条用颜色区分

cnetplot(erich.go.BP,showCategory = 5,circular=T,  ###画为圈图
         colorEdge=T)

library(ggplot2)
library(dplyr)

plt_DEG <- function(result){
  DEG <- data.frame()
  for (i in 1:9){
    gene <- unlist(strsplit(result$geneID[i],"/"))
    ti <- result$Description[i]

    DEG_go_neg <- DEG_matrix %>% dplyr::select(mgi_symbol,logFC) %>% filter(mgi_symbol %in% gene & logFC<0 ) %>%  mutate(type = "down",group = ti)
    DEG_go_pos <- DEG_matrix %>% dplyr::select(mgi_symbol,logFC)  %>% filter(mgi_symbol %in% gene & logFC>0 ) %>% mutate(type = "up", group = ti)
    DEG_go <- rbind(DEG_go_neg,DEG_go_pos)
    DEG <-rbind(DEG, DEG_go)

  }
  ggplot(DEG, aes(x = mgi_symbol, y = logFC, color = type)) + 
    geom_point() + 
    coord_flip() +
    facet_wrap( ~ group, ncol = 3,scales = "free_y") + 
    theme_grey() + 
    ggtitle("Differentially Expressed Gene") +
    theme(plot.title = element_text(hjust = 0.5))
}

plt_DEG(erich.go.BP@result)