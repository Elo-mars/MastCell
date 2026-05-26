Seurat_reso = FC2
head(Seurat_reso@meta.data)
min(Seurat_reso@meta.data$nCount_RNA)
min(Seurat_reso@meta.data$nFeature_RNA) 
max(Seurat_reso@meta.data$percent.mt) 

# [1] 1000.008
# [1] 401
# [1] 19.99997


RESO="RNA_snn_res.1.4"
u=Seurat_reso@meta.data
write.table(u, "Integrated_metadata.txt",col.names=NA, sep="\t")
Seurat_reso <- SetIdent(object = Seurat_reso, value = RESO)
table(Seurat_reso@active.ident)

plots <- FeaturePlot(Seurat_reso, features="nFeature_RNA", cols = c("cadetblue2", "darkred"), pt.size=0.1)
plots + scale_colour_continuous("# genes", low="cadetblue2", high="darkred") + ggtitle("Number of genes detected") + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15), plot.title = element_text(size=15))
ggsave("graph_nFeature.pdf", height = 5, width = 5)

plots <- FeaturePlot(Seurat_reso, features="nCount_RNA", cols = c("cadetblue2", "darkred"), pt.size=0.1)
plots + scale_colour_continuous("# UMIs", low="cadetblue2", high="darkred") + ggtitle("Number of UMIs detected") + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15), plot.title = element_text(size=15))
ggsave("graph_nCount.pdf", height = 5, width = 5)

plots <- FeaturePlot(Seurat_reso, features="percent.mt", cols = c("cadetblue2", "darkred"), pt.size=0.1)
plots + scale_colour_continuous("% mito genes \ncontent", low="cadetblue2", high="darkred") + ggtitle("% mito gene content") + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15), plot.title = element_text(size=15))
ggsave("graph_percentmito.pdf", height = 5, width = 5)

plots <- FeaturePlot(Seurat_reso, features="percent.ribo", cols = c("cadetblue2", "darkred"), pt.size=0.1)
plots + scale_colour_continuous("% ribo genes \ncontent", low="cadetblue2", high="darkred") + ggtitle("% ribo gene content") + theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15), plot.title = element_text(size=15))
ggsave("graph_percentribo.pdf", height = 5, width = 5)

head(Seurat_reso@meta.data)
table(Seurat_reso@active.ident)

features=FetchData(Seurat_reso, vars = c("nFeature_RNA","nCount_RNA", "percent.mt","percent.ribo", RESO), cells = NULL)
write.table(features, "table_features.txt",col.names=NA, sep="\t")

f=ordered(features$RNA_snn_res.1.4 , levels=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38))

pdf("boxplot_nFeature.pdf")
box = boxplot(nFeature_RNA ~ f , data=features, xlab="# genes", ylab="Cluster", main="Number of Genes per cluster" ,col = "dodgerblue3",  horizontal=TRUE)
dev.off()

pdf("boxplot_nCount.pdf")
box = boxplot(nCount_RNA ~ f , data=features, xlab="# UMIs", ylab="Cluster", main="Number of UMIs per cluster" ,col = "dodgerblue3",  horizontal=TRUE)
dev.off()

pdf("boxplot_percentmito.pdf")
box = boxplot(percent.mt ~ f , data=features, xlab="% mito", ylab="Cluster", main="Percentage of mito genes per cluster" ,col = "dodgerblue3",  horizontal=TRUE)
dev.off()

pdf("boxplot_percentribo.pdf")
box = boxplot(percent.ribo ~ f , data=features, xlab="% ribo", ylab="Cluster", main="Percentage of ribo genes per cluster" ,col = "dodgerblue3",  horizontal=TRUE)
dev.off()


FC=FC2

Idents(FC)=RESO
table(FC@active.ident)

# # check visually in which clusters you have high % mito and low %ribo
# Clusters: 19,11,27,29 (parts)
# ranges:
# sub= subset(FC, idents = c("19","11","27","29")) #  cells
#  head(sub@meta.data,2)
#  median(sub@meta.data$percent.mt)  
#  [1]  14.10955
#  median(sub@meta.data$percent.ribo)  
#  [1] 6.780482
#  
#  try: 
# keep: pct mito <14.1 and  ribo >6.8 
FC
Idents(FC) = "RNA_snn_res.1.4"
reso="1.4"
clusters="39"

pdf(paste0("Filter_dead_",reso,"_",clusters,"clusters_PRE.pdf"))
UMAPPlot(FC, reduction = "umap", pt.size=0.1,label=T, label.size = 5) + labs(title = paste0("UMAP with res ",reso))+ theme(plot.title = element_text(hjust = 0.5))
dev.off()


FCfiltered <- subset(FC, subset = percent.mt < 14.1 & percent.ribo > 6.8) 
FCfiltered

Idents(FCfiltered) = "RNA_snn_res.1.4"
reso="1.4"
clusters="39"

pdf(paste0("Filter_dead_",reso,"_",clusters,"clusters_POST.pdf"))
UMAPPlot(FCfiltered, reduction = "umap", pt.size=0.1,label=T, label.size = 5) + labs(title = paste0("UMAP with res ",reso))+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

ORGANISM="Human"
FC=FCfiltered

saveRDS(FC,"rds_nodead_noDoublets_292315cells.rds")
#for (resolution in c("0.4","0.6","0.8","1","1.2","1.4","1.6")){
resolution="1.4"
Seurat_reso <- SetIdent(object = FC, value = paste0('RNA_snn_res.',resolution))
library(future)
plan("multicore", workers = 8)   
options(future.globals.maxSize = 100000 * 1024^2) 

if(ORGANISM=="Mouse"){                                       
    genes.use <- grep(pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",rownames(Seurat_reso), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
} else { 
    genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",rownames(Seurat_reso), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
}


# ======================================================== #
# FindAllMarkers - generate gene lists for each resolution #
# ======================================================== #

# Make sure you are in an environment using Seurat 4! not 5, the algorithm for calculating FAM changed a lot between these 2 versions and the results with Seurat 4 are preferred by us
ORGANISM="Human"
for (resolution in c("0.4","0.6","0.8","1","1.2","1.4","1.6")){
Seurat_reso <- SetIdent(object = FC, value = paste0('RNA_snn_res.',resolution))
library(future)
plan("multicore", workers = 60)   
options(future.globals.maxSize = 100000 * 1024^2) 

if(ORGANISM=="Mouse"){                                       
    genes.use <- grep(pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",rownames(Seurat_reso), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
} else { 
    genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",rownames(Seurat_reso), value=TRUE, invert=TRUE) #get list of non-ribosomal genes 
}


Seurat.markers <- FindAllMarkers(Seurat_reso,test.use="wilcox",  min.pct = 0.1, logfc.threshold = 0.25, features=genes.use) # no ribo genes
write.table(Seurat.markers, paste0("Combined_res.",resolution,"_wilcox_markers_unannotated_Clusters.txt"),col.names=NA, sep="\t")
tsne.markers=Seurat.markers         
topn = tsne.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
topn$gene                                                     
pdf(paste0("Combined_res.",resolution,"_heatmap_top10genes.pdf"),width=25, height=20)                                                                                                                                                                                     
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 10 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
topn$gene
pdf(paste0("Combined_res.",resolution,"_heatmap_top15genes.pdf"),width=25, height=25)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 15 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
topn$gene                                                                                                                             
pdf(paste0("Combined_res.",resolution,"_heatmap_top20genes.pdf"),width=25, height=30)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 20 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)  
topn$gene
pdf(paste0("Combined_res.",resolution,"_heatmap_top5genes.pdf"),width=25, height=15)
print(DoHeatmap(object = Seurat_reso, features=topn$gene, group.by = "ident", label = TRUE, angle=45) + ggtitle("top 5 genes expressed per cluster") + NoLegend())
dev.off()
topn = tsne.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC) 
topn$gene
topn
write_csv(topn, paste0("Combined_res.",resolution,"_top20genespercluster_logFC.csv"))
}

