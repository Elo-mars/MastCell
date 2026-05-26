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

# ============== #
# Doublet Finder #
# ============== #

# DOUBLET FINDER 3 
doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE) {
  require(Seurat); require(fields); require(KernSmooth)
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }
    
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
}

Object = combined
Object = SetIdent(object = Object, value = "orig.ident")
table(Object@active.ident)
head(Object@meta.data)
Object$pANN <- "NA"
Object$pANNPredictions <- "NA"
head(Object@meta.data)

for(sample in unique(Object$orig.ident)){
  sample.cluster <- subset(Object, idents = sample)
  print(paste0("sample:", sample))
  length(rownames(sample.cluster@meta.data))
  expected.doublets <- ceiling(0.039 * length(rownames(sample.cluster@meta.data)))
  sample.cluster  <- doubletFinder_v3(sample.cluster, PCs = 1:20, pN = 0.25, pK = 0.01, nExp = expected.doublets, reuse.pANN = FALSE, sct=TRUE)
  sample.cluster$pANN <- sample.cluster@meta.data[colnames(sample.cluster), paste("pANN_0.25_0.01", expected.doublets, sep = "_")]
  sample.cluster$pANNPredictions <- sample.cluster@meta.data[colnames(sample.cluster), paste("DF.classifications_0.25_0.01", expected.doublets, sep = "_")]
  Object$pANN[colnames(sample.cluster)] <- sample.cluster$pANN[colnames(sample.cluster)]
  Object$pANNPredictions[colnames(sample.cluster)] <- sample.cluster$pANNPredictions[colnames(sample.cluster)]
  sample.cluster <- NULL
}

head(Object@meta.data)
combined@meta.data = Object@meta.data
head(combined@meta.data)

x = combined@meta.data
write.table(x, "DF_metadata.txt",col.names=NA, sep="\t")

pdf("doublets.pdf")
DimPlot(Object,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "pANNPredictions" )+theme(aspect.ratio = 1)
dev.off()

Idents(combined)="pANNPredictions"
table(combined@active.ident)

#subset all the doublets
combined
FC2 = subset(combined, idents="Singlet")
FC2 

# umap only these
pdf("singlets_per_PANN.pdf")
DimPlot(FC2,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "pANNPredictions" )+theme(aspect.ratio = 1)
dev.off()

pdf("singlets_per_cluster.pdf")
DimPlot(FC2,pt.size = 0.1,label=F, label.size = 0,reduction = "umap",group.by = "RNA_snn_res.1.4" )+theme(aspect.ratio = 1)
dev.off()

saveRDS(FC2, "FC_removed_doublet_cells.rds")
combined=FC2

p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident", pt.size=0.1, label = F, label.size=3.5)
p2 <- DimPlot(combined, reduction = "umap", group.by = "Status",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 
p4 <- DimPlot(combined, reduction = "umap", group.by = "Subtype",pt.size=0.1, label = F, label.size=3.5,  repel = TRUE) 

pdf("umap_combined.pdf", width=15, height=8)
plot_grid(p1, p2,  p4)
dev.off()

saveRDS(combined, file = "combined_FC.rds")

