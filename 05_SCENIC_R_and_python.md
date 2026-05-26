
```r
library(Seurat)
PROJECT="/path/Figure_IBS"
setwd(PROJECT)
seu_object = readRDS("Figure_IBS.rds")

# if you need to downsample
Idents(seu_object) <- 'Cluster'
cell.list <- WhichCells(seu_object, idents = levels(seu_object@meta.data$Cluster), downsample = 4000) # here, for example, we go to 4000 cells per cluster max

downsampled <- seu_object[, cell.list]

dir.create("SCENIC")
setwd(paste0(PROJECT,"/SCENIC"))
saveRDS(downsampled, "downsampled.rds") # you'll need it later
seu_object=downsampled

exprMat <- BiocGenerics::t(GetAssayData(object = seu_object[["RNA"]])) 
head(exprMat)
x = data.frame("cell_id"=rownames(exprMat),exprMat)
x[0:5,0:5] # sanity check
write.table(data.frame("cell_id"=rownames(exprMat),exprMat),file="exprMat.tsv",row.names=FALSE,col.names=TRUE,sep="\t",quote = FALSE)

q()
n
```

SCENIC: check there for info: https://pyscenic.readthedocs.io/en/latest/installation.html

```python
SCENIC_IMAGE=/path/software/SCENIC/aertslab-pyscenic-0.12.1.sif
SCENIC_RESOURCES=/path/software/SCENIC
DATA_FOLDER=$1
NR_WORKERS=25

singularity exec -B $DATA_FOLDER:/scenicdata,$SCENIC_RESOURCES:/resources \
$SCENIC_IMAGE \
pyscenic grn \
--num_workers $NR_WORKERS \
-o /scenicdata/expr_mat.adjacencies.tsv \
/scenicdata/exprMat.tsv \
/resources/human_tfs.txt


singularity exec -B $DATA_FOLDER:/scenicdata,$SCENIC_RESOURCES:/resources \
$SCENIC_IMAGE \
pyscenic ctx \
/scenicdata/expr_mat.adjacencies.tsv \
/resources/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather \
/resources/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /scenicdata/exprMat.tsv \
--mode "dask_multiprocessing" \
--output /scenicdata/regulons.csv \
--num_workers $NR_WORKERS


singularity exec -B $DATA_FOLDER:/scenicdata,$SCENIC_RESOURCES:/resources \
$SCENIC_IMAGE \
pyscenic aucell \
/scenicdata/exprMat.tsv \
/scenicdata/regulons.csv \
-o /scenicdata/auc_mtx.csv \
--num_workers $NR_WORKERS
```

```r
library(Seurat)
PROJECT="/path/Figure_IBS/SCENIC/"
setwd(PROJECT)
seu_object = readRDS("downsampled.rds")
output_path=paste0(PROJECT)

scenic_results <- read.table(paste0(output_path,"auc_mtx.csv"), header=T, sep=",",row.names=1)
colnames(scenic_results) <- sub("...$","",colnames(scenic_results))
scenic_results <- t(scenic_results)
scenic_results <- scenic_results[,colnames(seu_object)]
seu_object[["SCENIC"]] <- CreateAssayObject(scenic_results)
seu_object <- ScaleData(seu_object, assay="SCENIC")
seu_object <- RunPCA(seu_object, features = rownames(seu_object), reduction.name = "pca_scenic", reduction.key = "pca_scenic_", assay="SCENIC", verbose = FALSE)
scenic.markers <- FindAllMarkers(seu_object, assay = "SCENIC",test.use = "roc", logfc.threshold = FALSE)
result <- list("seu"=seu_object, "markers"= scenic.markers)
write.table(scenic.markers, "SCENIC.txt",col.names=NA, sep="\t")

library(viridis)

PlotScenicMarkers <- function(seu_object, markers, output_path){
  for(c in unique(markers$cluster)){
    out <- paste0(output_path,c)
    dir.create(out)
    genes <- markers[markers$cluster == c,"gene"]
    for(gi in 1:length(genes)){
        g=genes[gi]
        DefaultAssay(seu_object) <- "SCENIC"
        pdf(file=paste0(out,"/",g,".pdf"),height=20,width=25)
        p1 <- DimPlot(seu_object, reduction="umap", pt.size=1)
        p2 <- FeaturePlot(seu_object, features = g, min.cutoff = "q7",reduction = "umap", pt.size=1) + scale_color_viridis(option="viridis")
        p3 <- VlnPlot(seu_object, features = g,pt.size=0.5) 
        p4 <- RidgePlot(seu_object, features = g)
        print(CombinePlots(list(p1,p2,p3,p4),ncol=2))
        dev.off()
    }
  }
}

plots = PlotScenicMarkers(seu_object, scenic.markers, output_path)

PlotScenicHeatmap <- function(seu_object, markers, output_path){
  seu.small <- subset(seu_object, downsample = 100)
  pdf(file=paste0(output_path,"SCENIC_heatmap.pdf"),height=30,width=12)
  print(DoHeatmap(seu.small, features = unique(markers$gene), assay = "SCENIC", angle = 90) + NoLegend())
  dev.off()
}

plots2 = PlotScenicHeatmap(seu_object, scenic.markers, output_path)

# =================== #
# Dotplot of top 5 TF #
# =================== #

Seurat_Dotplot = Seurat_Final

DefaultAssay(Seurat_Dotplot) <- "SCENIC"
Seurat_Dotplot=SetIdent(Seurat_Dotplot, value="Cluster")

tsne.markers=read.table("SCENIC.txt", sep="\t", header=TRUE)
topn = tsne.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
topn$gene
topn
df = as.data.frame(topn)
df
Markers = df[,"gene"]
Markers = unique(Markers)
head(Markers)

scale="cividis"
markers.to.plot <- Markers
Idents(Seurat_Dotplot) = "Cluster"

plot = DotPlot(Seurat_Dotplot,assay="SCENIC", features = rev(markers.to.plot), dot.scale = 4) + scale_color_viridis(option= scale) + RotatedAxis() + theme(plot.margin = unit(c(0.1,0.1,0.1,2), "cm")) + theme(axis.text=element_text(size=10))
ggsave("TF_dotplot_top5_genes_all_clusters.pdf",height = 5, width = 15)
```

