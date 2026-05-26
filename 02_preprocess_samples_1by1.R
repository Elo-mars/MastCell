PROJECT="/path/project"
# sample name of the directory where the data is
SAMPLE="sample_name" 
samplename=SAMPLE
ORGANISM="Human"
name="PIN_xxx"
subDir<-paste0("Result_",SAMPLE)
dir.create(file.path(PROJECT,"/sample1_by_1/", subDir), showWarnings = FALSE)
opath_s<-paste0(PROJECT,"/sample1_by_1/",subDir,"/")
opath<-opath_s
setwd(opath)

library("Seurat")
library("plyr")
library("dplyr")
library("future")
library("ggplot2")
library("cowplot")
library("tidyverse")
library("grid")
library("gridExtra")
library("gtable")
library("readr")
library("beanplot")
library("gplots") 
library("RColorBrewer")
library("devtools")
library("sctransform")
library("harmony")
library("fields")
library("viridis")

##########################
#######  SoupX     #######
##########################                                                                                           
testDir<-paste0(PROJECT,"/samples/",SAMPLE,"/outs")
sc = load10X(testDir)
sc = autoEstCont(sc)
out = adjustCounts(sc)
sample= CreateSeuratObject(out)
saveRDS(sample, file =paste0(opath,name,"_After_SoupX.rds"))
### remove the genes expressed in less than 3 cells
gene_sum<-Matrix::rowSums(sample@assays$RNA@data>=3)
gene_keep_index<-(which(gene_sum>0))
count.data <- GetAssayData(object = sample[["RNA"]])[rownames(sample)[gene_keep_index],]
sample <- SetAssayData(object = subset(sample,features=rownames(sample)[gene_keep_index]), new.data = count.data,assay = "RNA")
saveRDS(sample, file =paste0(opath,name,"_After_SoupX_3cellgenes.rds"))

head(sample@meta.data)
sample@meta.data$orig.ident = name
sample@meta.data$Status <- "IBS"
sample@meta.data$Subtype <- "IBS-C"

# ============================================== #
# 1) QC and selecting cells for further analysis #
# ============================================== #

if(ORGANISM=="Mouse"){
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
} else {
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
}

# check ribosomal content
if(ORGANISM=="Mouse"){
sample<- PercentageFeatureSet(sample, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa", col.name = "percent.ribo", assay = 'RNA')
} else {
	sample<- PercentageFeatureSet(sample, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo", assay = 'RNA')
}

pdf(paste0(name,"_VlnPlot.pdf"), width = 15, height = 15)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0.5,ncol = 2)
dev.off()

# ====== #
# Filter #
# ====== #

y = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > -Inf, nFeature_RNA < +Inf, percent.mt > -Inf, percent.mt < +Inf, nCount_RNA > -Inf, nCount_RNA < +Inf) %>%
  summarise(n=n())
write.table(y, paste0(name,"_beforeQC_numberofcells.txt"),col.names=NA, sep="\t")

x = sample@meta.data %>% group_by(orig.ident) %>%
  filter(nFeature_RNA > 300, nFeature_RNA < 6000, percent.mt > -Inf, percent.mt < 30, nCount_RNA > 1000, nCount_RNA < 50000) %>%
  summarise(n=n())
 x
write.table(x, paste0(name,"_QC_numberofcells.txt"),col.names=NA, sep="\t")

# features are genes and count_RNA are UMIs
sample <- subset(sample, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 30 & nCount_RNA > 1000 & nCount_RNA < 50000)


fileConn<-file("FILTERS.txt")
writeLines(c("nfeatures","300-6.000","numis","1.000-50.000","%mt","<30"), fileConn)
close(fileConn)


pdf(paste0(name,"_filtered_VlnPlot.pdf"), width = 15, height = 8)
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
dev.off()

# ==================== #
# Normalizing the data #
# ==================== #

# After removing unwanted cells from the dataset, the next step is to normalize the data.
Norm_sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)

# ============================================================== #
# Identification of highly variable features (feature selection) #
# ============================================================== #

FHF <- FindVariableFeatures(Norm_sample, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(FHF), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(FHF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(paste0(name,"_Highly_Variable_genes.pdf"), width=10, height=5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

all.genes=rownames(FHF)

library(future)
plan("multiprocess", workers = 60)
options(future.globals.maxSize = 10000 * 1024^2) 

# ===== #
# Scale #
# ===== #

Scaled <- ScaleData(FHF, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mt"))

# ======= #
# Run PCA #
# ======= #

library(future)
plan("multiprocess", workers = 60) 
options(future.globals.maxSize = 10000 * 1024^2) 
Scaled <- RunPCA(Scaled, features = VariableFeatures(Scaled), npcs = 30, ndims.print = 1:5, nfeatures.print = 30)

pdf(paste0(name,"_PCA.pdf"))
PCAPlot(object = Scaled, dims=c(1, 2)) 
dev.off()

# ======================================================== #
# Determine statistically significant principal components #
# ======================================================== #

Scaled <- JackStraw(Scaled, num.replicate = 100, dims = 30)
Scaled <- ScoreJackStraw(Scaled, dims = 1:30)

pdf(paste0(name,"_JackStraw.pdf"))
JackStrawPlot(Scaled, dims = 1:30)
dev.off()
pdf(paste0(name,"_Elbow.pdf"))
ElbowPlot(Scaled, ndims = 30, reduction = "pca")
dev.off()

NPC=30

# ==================== #
# 3) Cluster the cells #
# ==================== #

FN <- FindNeighbors(Scaled, reduction = "pca", dims = 1:NPC, k.param = 30, compute.SNN = TRUE, nn.method = "rann", nn.eps = 0)
plan("multiprocess", workers = 60) 
options(future.globals.maxSize = 10000 * 1024^2) 
FC <- FindClusters(FN, resolution = c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.4,1.6))

saveRDS(FC,paste0(samplename,"_FC.rds"))
head(FC@meta.data)

# ================================================ #
# Run non-linear dimensional reduction (UMAP/tSNE) #
# ================================================ #

Idents(FC)

# Color by sample
Seurat_per_sample_Ident <- SetIdent(object = FC, value = 'orig.ident')
Seurat_per_sample <- RunUMAP(Seurat_per_sample_Ident, dims = 1:NPC)
Seurat_per_sample

pdf(paste0(name,"_Plot_per_sample.pdf"))
UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = "UMAP per sample")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(name,"_Plot_per_sample_nolegend.pdf"))
UMAPPlot(Seurat_per_sample, reduction = "umap", pt.size=0.1,label=FALSE) +labs(title = "UMAP per sample")+ theme(plot.title = element_text(hjust = 0.5)) +NoLegend()
dev.off()

head(FC@meta.data)

pdf(paste0(name,"_Clusters_diff_res.pdf"))
Seurat <- SetIdent(object = FC, value = 'RNA_snn_res.1')
Seurat <- RunUMAP(Seurat, dims = 1:NPC)
UMAPPlot(Seurat, reduction = "umap", pt.size=0.25,label=T) +labs(title = "UMAP with res 0.4")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()
