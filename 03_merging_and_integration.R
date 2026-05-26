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
library("harmony")                                                                                                                                                                                                                                                              
library("fields")                                                                                                                                                                                                                                                               
sessionInfo()                                                                                                                                                                                                                                                                   

PROJECT="/path/here" 

setwd(PROJECT)

SAMPLElist = c("PIN_xxxx","PIN_yyyy","PIN_zzzz",...)                                                                                            

sample_list<-list()                                                                                                                  
for (i in (1:length(SAMPLElist))) {                                                                                                                                                                                                                                         
  sample=SAMPLElist[i]                                                                                                                                                                                                                                                      
  if (file.exists(paste0(PROJECT,"/",sample,".rds"))){                                                                               
    so<-readRDS(paste0(PROJECT,"/",sample,".rds"))                                                                                                                                                                                                                          
      sample_list[[i]]<-so                                                                                                                                                                                                                                                  
}}                                                                                                                                   


GR <- merge(x = sample_list[[1]],y=c(sample_list[[2]],sample_list[[3]],sample_list[[4]],sample_list[[5]],sample_list[[6]],sample_list[[7]],sample_list[[8]],sample_list[[9]],sample_list[[10]],sample_list[[11]],sample_list[[12]],sample_list[[13]],sample_list[[14]],sample_list[[15]],sample_list[[16]],sample_list[[17]],sample_list[[18]],sample_list[[19]],sample_list[[20]],sample_list[[21]],sample_list[[22]],sample_list[[23]],sample_list[[24]]),add.cell.ids = SAMPLElist, project = PROJECT)
x = GR@meta.data %>% group_by(orig.ident) %>%filter(nFeature_RNA > 400, nFeature_RNA < 6000, percent.mt > -Inf, percent.mt < 20, nCount_RNA > 1000, nCount_RNA < 50000) %>% summarise(n=n())                                                
x                                                                
write.table(x, "GR_QC_numberofcells.txt",col.names=NA, sep="\t")

GR <- subset(GR, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20 & nCount_RNA > 1000 & nCount_RNA < 50000) 

fileConn<-file("FILTERS.txt")
writeLines(c("nfeatures","400-6.000","numis","1.000-50.000","%mt","<20"), fileConn)
close(fileConn)

# ====== #
# Violin #
# ====== #
Idents(GR) = "orig.ident"
pdf("ALL_VlnPlot.pdf",width = 15, height = 15)
VlnPlot(GR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0.5,ncol = 2)
dev.off()

GR2=GR
Idents(GR2) = "GC"
pdf("ALL_VlnPlot_GC.pdf",width = 25, height = 15)
VlnPlot(GR2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0.5,ncol = 2)
dev.off()

combined = GR
tail(combined@meta.data)
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

combined

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
print("FindVarFeatures done")

combined <- ScaleData(combined, do.scale = TRUE,do.center=TRUE,  features = VariableFeatures(object = combined), vars.to.regress = c("nCount_RNA", "percent.mt"))
print("Scale done") 

combined = RunPCA(combined,features = VariableFeatures(object = combined), npcs = 35, ndims.print = 1:5, nfeatures.print = 30)
print("RunPCA done")
saveRDS(combined, "before_runharmony.rds")
combined <- RunHarmony(combined, c("orig.ident"))
saveRDS(combined, "after_runharmony.rds")

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:35)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:35, k.param = 30, compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) 
combined <- FindClusters(combined, resolution = c(0.4,0.6,0.8,1,1.2,1.4))

head(combined@meta.data)

p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident", pt.size=0.1, label = F, label.size=3.5)
p2 <- DimPlot(combined, reduction = "umap", group.by = "Status",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 
p4 <- DimPlot(combined, reduction = "umap", group.by = "Subtype",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 

pdf("umap_combined.pdf", width=15, height=8)
plot_grid(p1, p2,  p4)
dev.off()

saveRDS(combined, file = "combined_FC.rds")

