library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(scater)

PROJECT="date_CD45_HV_vs_IBS-D"

mainDir<-"/path/milo"
subDir<-paste0("miloR_",PROJECT)
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
opath<-paste0(mainDir,"/",subDir,"/")
rdspath<-paste0(mainDir,"/")
setwd(opath)

# Really needs to be a Seurat4 object
FC<-readRDS(file=paste0(rdspath,'obj.rds')) 
condition="HV_vs_IBSD"
print(condition)
Idents(FC)<-"Subtype"
table(FC$Subtype)
FC$sample<-FC$orig.ident

FC@meta.data$Subtype[FC@meta.data$Subtype == "HV"] <- "2_HV"
FC@meta.data$Subtype[FC@meta.data$Subtype == "IBS-D"] <- "1_IBS-D"
FC@meta.data$Subtype[FC@meta.data$Subtype == "IBS-C"] <- "1_IBS-C"
FC@meta.data$Subtype[FC@meta.data$Subtype == "IBS-M"] = "1_IBS-M"
Idents(FC)<-"Subtype"
table(FC@meta.data$orig.ident, FC@meta.data$Subtype)
sub_SO<-FC[,!(FC@active.ident %in% c("1_IBS-C","1_IBS-M"))] # the one you remove

table(sub_SO$Subtype)
table(sub_SO$orig.ident)
sub_SO_se <- as.SingleCellExperiment(sub_SO)

sub_SO_se$sample<-sub_SO_se$orig.ident
# DA testing
sub_SO_milo <- Milo(sub_SO_se)
sub_SO_milo <- buildGraph(sub_SO_milo, k = 35, d = 35, reduced.dim = "PCA")
sub_SO_milo <- makeNhoods(sub_SO_milo, prop = 0.1, k = 35, d=35, refined = TRUE, reduced_dims = "PCA")
sub_SO_milo <- countCells(sub_SO_milo, meta.data = as.data.frame(colData(sub_SO_milo)), sample="sample")
head(nhoodCounts(sub_SO_milo))
sub_SO_milo_design <- data.frame(colData(sub_SO_milo))[,c("sample", "Subtype")]

## Convert batch info from integer to factor
sub_SO_milo_design$Subtype <- as.factor(sub_SO_milo_design$Subtype) 
dim(sub_SO_milo_design)
sub_SO_milo_design <- distinct(sub_SO_milo_design)
dim(sub_SO_milo_design)
rownames(sub_SO_milo_design) <- sub_SO_milo_design$sample

sub_SO_milo_design 

plan("multisession", workers = 10)
options(future.globals.maxSize = 500000 * 1024^2)

sub_SO_milo <- calcNhoodDistance(sub_SO_milo, d=35, reduced.dim = "PCA")
sub_SO_milo 

da_results <- testNhoods(sub_SO_milo, design = ~ Subtype, design.df = sub_SO_milo_design)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# Inspecting DA testing results
pdf(paste0(condition,"_pvalue_histogram.pdf"))
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
dev.off()

pdf(paste0(condition,"_pvalue_hline.pdf"))
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)
dev.off()

sub_SO_milo <- buildNhoodGraph(sub_SO_milo)
umap_pl<-DimPlot(sub_SO, split.by='Subtype',group.by='Subtype',cols=c("#a31818", "#2171b5"),ncol=1)

# Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(sub_SO_milo, da_results, layout="UMAP",alpha=0.1, pt.size=1) 

pdf(paste0(condition,"_umap_nl_graph.pdf"),width=6, height=10)
umap_pl + nh_graph_pl + plot_layout(guides="collect")
dev.off()

# color groups
da_results <- groupNhoods(sub_SO_milo, da_results, max.lfc.delta = 10, overlap = 1) 
head(da_results)

pdf(paste0(condition,"_plotNhoodGroups.pdf"),height= 8, width=9)
plotNhoodGroups(sub_SO_milo, da_results, show_groups = NULL, alpha=0.05, delta=3)
dev.off()

da_results <- annotateNhoods(sub_SO_milo, da_results, coldata_col = "Cluster")
head(da_results)

pdf(paste0(condition,"_Cluster_fraction.pdf"))
ggplot(da_results, aes(Cluster_fraction)) + geom_histogram(bins=50)
dev.off()

da_results$anno_celltype_cluster <- ifelse(da_results$Cluster_fraction < 0.7, "Mixed", da_results$Cluster)

pdf(paste0(condition,"_plotDAbeeswarm.pdf"),height= 15, width=17)
plotDAbeeswarm(da_results, group.by = "Cluster") + theme(text = element_text(size = 40)) + scale_colour_gradient2(limits = c(-10, 9)) 
dev.off()
write.table(da_results, paste0(condition,"_da_results.txt"),col.names=NA, sep="\t")

library(dplyr)
MILO_results=da_results
head(MILO_results)
z <- MILO_results %>% group_by(anno_celltype_cluster) %>% summarise(across(everything(), mean))
z
write.table(z, paste0("Milo_",PROJECT,"_",condition,"_mean.txt"),col.names=NA, sep="\t")

