# create an environment only for this
source activate Xenium
R 

PROJECT="/path/Xenium"
NEW <- "date_Xenium"
main_folder_path <- file.path(PROJECT,NEW)
main_folder_path
dir.create(main_folder_path)
setwd(main_folder_path)
getwd()

library(ggplot2)
library(dplyr)
library(sp)
library(data.table)
library(Seurat)
library(patchwork)
library(magrittr)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize= 8912896000)

plan("multicore", workers = 15) 
options(future.globals.maxSize = 400000 * 1024^2) 

 load_xenium<-function(folder, name){
     data <- ReadXenium(data.dir =paste0(folder, '/', sep=""), 
                        type = c("centroids", "segmentations"))
    assay <- "Xenium"
   segmentations.data <- list(
     "centroids" = CreateCentroids(data$centroids),
     "segmentation" = CreateSegmentation(data$segmentations)
   )
   coords <- CreateFOV(
     coords = segmentations.data,
     type = c("segmentation", "centroids"),
     molecules = data$microns,
     assay = assay
   )
   
   xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
   xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
   xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
   xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
   fov <- name
   
   pdf(paste0("Before_QC_Vln_",name,'.pdf'),width=12, height=7)
   print(VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0.1))
   dev.off()
   
   xenium.obj[[fov]] <- coords
   xenium.obj<-subset(xenium.obj, subset=nFeature_Xenium >1 & nCount_Xenium >10) 
   xenium.obj$technique<-"Xenium"
   
   pdf(paste0("After_QC_Vln_",name,'.pdf'),width=12, height=7)
   print(VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0.1))
   dev.off()

   xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", clip.range = c(-10, 10))
   xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
   xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
   xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
   xenium.obj <- FindClusters(xenium.obj, resolution = 0.4)
   xenium.obj$slide<-name
   print(head(xenium.obj@meta.data,3))
   saveRDS(xenium.obj, file = paste0(name,'.rds', sep=""))
 }
 
 dirs<-list.dirs('/staging/leuven/stg_00075/Project/240806_Hind_Xenium/SLIDES', recursive = F)
 names<-c("PIN1561","PIN1635","PIN1825b","PIN1861","PIN1913","PIN1914","PIN1916","PIN4063")
 for(i in 1:length(dirs)){
   load_xenium(folder = dirs[i], name = names[i])
 }

 setwd(main_folder_path)
 PIN1561=readRDS("/.../PIN1561.rds")
 PIN1825b=readRDS("/.../PIN1825b.rds")
 PIN4063=readRDS("/.../PIN4063.rds")
 PIN1913=readRDS("/.../PIN1913.rds")
 PIN1914=readRDS("/.../PIN1914.rds")
 PIN1635=readRDS("/.../PIN1635.rds")
 PIN1916=readRDS("/.../PIN1916.rds")
 PIN1861=readRDS("/.../PIN1861.rds")
 
  head(PIN1916@meta.data,2)
 
 #               orig.ident nCount_Xenium nFeature_Xenium nCount_BlankCodeword
 # aaaallln-1 SeuratProject            21              16                    0
 # aaabpjgh-1 SeuratProject            23              21                    0
 #            nFeature_BlankCodeword nCount_ControlCodeword
 # aaaallln-1                      0                      0
 # aaabpjgh-1                      0                      0
 #            nFeature_ControlCodeword nCount_ControlProbe nFeature_ControlProbe
 # aaaallln-1                        0                   0                     0
 # aaabpjgh-1                        0                   1                     1
 #            technique nCount_SCT nFeature_SCT SCT_snn_res.0.4 seurat_clusters
 # aaaallln-1    Xenium         63           18              14              14
 # aaabpjgh-1    Xenium         59           23               8               8
 #              slide
 # aaaallln-1 PIN1916
 # aaabpjgh-1 PIN1916
 
 names<-c("PIN1561","PIN1635","PIN1825b","PIN1861","PIN1913","PIN1914","PIN1916","PIN4063")
 Xenium.big <- merge(PIN1561, y =c(PIN1635,PIN1825b,PIN1861,PIN1913,PIN1914,PIN1916,PIN4063), add.cell.ids =names, project = "PROJECT")
 DefaultAssay(Xenium.big) <- "Xenium"
 Xenium.big
  
 plan("multicore", workers = 15) 
 options(future.globals.maxSize = 400000 * 1024^2) 
 
 Xenium.big <-SCTransform(Xenium.big, assay = "Xenium", variable.features.n=422) 
 Xenium.big <- RunPCA(Xenium.big, npcs = 30)
 Xenium.big <- RunUMAP(Xenium.big, dims = 1:30)
 Xenium.big <- FindNeighbors(Xenium.big, reduction = "pca", dims=1:30)
 Xenium.big <- FindClusters(Xenium.big, resolution = 0.4)
 head(Xenium.big@meta.data,3)
 
 
 pdf(paste0("merged_DimPlot_slide.pdf"), width=20, height=20)
 DimPlot(Xenium.big, label = T, label.size = 3, reduction = "umap", group.by = "slide", shuffle = T) #+NoLegend()
 dev.off()
 
 pdf(paste0("merged_DimPlot_automaticclusters.pdf"),width=20, height=20)
 DimPlot(Xenium.big, label = T, label.size = 15, reduction = "umap", group.by = "seurat_clusters", shuffle = T) #+NoLegend()
 dev.off()
 
 Xenium.big@meta.data$orig.ident=""
 Xenium.big@meta.data$orig.ident=Xenium.big@meta.data$slide
 head(Xenium.big@meta.data,2)
 tail(Xenium.big@meta.data,2)
 
 names<-c("PIN1561","PIN1635","PIN1825b","PIN1861","PIN1913","PIN1914","PIN1916","PIN4063")
 for (i in names){
   DefaultBoundary(Xenium.big[[i]]) <- "segmentation" #prefered over centroids
 }
 
 saveRDS(Xenium.big,"MERGED.rds")

# we will annotate thanks to our SC object
SC=readRDS("/.../FINAL.rds")
head(SC@meta.data,2)
u=table(SC@active.ident)
u
write.table(u, "cells_per_cluster_SConly.txt",col.names=NA, sep="\t")

pdf("SC.pdf")
DimPlot(SC, pt.size=2)
dev.off()


# from SC, take your count matrix, subset with ony genes that are also found in your Xenium
# extract too the metadata
# with the new count matrix and the metadata, create a new SC object which is also only with 422 genes
DefaultAssay(SC)<-"RNA"
seu <- CreateSeuratObject(counts = SC@assays$RNA$counts, assay = "Xenium")
seu@meta.data<-SC@meta.data
head(seu@meta.data,2)

length(rownames(Xenium.big)) # 422 genes
length(rownames(seu)) # 17026 genes
head(rownames(Xenium.big))
head(rownames(seu))

spatial_features <- intersect(rownames(Xenium.big), rownames(seu))
head(spatial_features)
length(spatial_features) # 413 genes

write.table(spatial_features, "spatial_features.txt",col.names=NA, sep="\t")

# look if there is duplicates
unique(spatial_features)

subsetted_seu <- subset(seu, features = spatial_features) # same cell # but with only 413 genes
subsetted_seu<-JoinLayers(subsetted_seu, assay  = "Xenium")
subsetted_seu<-NormalizeData(object = subsetted_seu, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Xenium")
subsetted_seu$technique<-"scRNA"
head(subsetted_seu@meta.data,2)
subsetted_seu=SCTransform(subsetted_seu, assay = "Xenium", variable.features.n =421)
subsetted_seu<- RunPCA(subsetted_seu, npcs = 30, features = rownames(subsetted_seu))
subsetted_seu <- RunUMAP(subsetted_seu, dims = 1:30)
DefaultAssay(subsetted_seu)<-"Xenium"
saveRDS(subsetted_seu, "reduced_SC.rds")
subsetted_seu

# same for Xenium, because not all the genes in your panel are in your single cell data, sometimes, they use of the shelf probes which are absolutely not in your SC
Xenium.big
DefaultAssay(Xenium.big)<-"SCT"
Xenium.big<-JoinLayers(Xenium.big, assay  = "Xenium")
Xen <- CreateSeuratObject(counts = Xenium.big@assays$Xenium$counts, assay = "Xenium")
Xen@meta.data<-Xenium.big@meta.data
head(Xen@meta.data,2)
Xen
  
downed_merged_Xenium <- subset(Xen, features = spatial_features) # still 716,380 cells with 413 genes
downed_merged_Xenium<-JoinLayers(downed_merged_Xenium, assay  = "Xenium")
downed_merged_Xenium<-NormalizeData(object = downed_merged_Xenium, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Xenium")
downed_merged_Xenium$technique<-"Xen"
head(downed_merged_Xenium@meta.data,2)
downed_merged_Xenium=SCTransform(downed_merged_Xenium, assay = "Xenium", variable.features.n =421)
downed_merged_Xenium<- RunPCA(downed_merged_Xenium, npcs = 30, features = rownames(downed_merged_Xenium))
downed_merged_Xenium <- RunUMAP(downed_merged_Xenium, dims = 1:30)
DefaultAssay(downed_merged_Xenium)<-"Xenium"
saveRDS(downed_merged_Xenium, "reduced_Xenium.rds")
downed_merged_Xenium


# ===================== #
# MY 2 OBJECTS TO MERGE #
# ===================== #

SC=subsetted_seu
SPATIAL=downed_merged_Xenium

# ==== #
# DotR #
# ==== #

# Single cell: CD45, HV only, reduced to the intersectrion of genes between single cell and Xenium
subsetted_seu
# Xenium: reduced to the intersectrion of genes between single cell and Xenium
downed_merged_Xenium

dim(subsetted_seu)
dim(downed_merged_Xenium)

plot_data_ref=as.data.frame(SC[["umap"]]@cell.embeddings)
head(plot_data_ref)

#                                                    umap_1     umap_2
# PIN_1194_GC114074_SI-GA-E10_AAAGAACCACGGCCAT-1 -3.8958291 -1.8868753
# PIN_1194_GC114074_SI-GA-E10_AAATGGATCAGTGATC-1 -3.1601912 -2.6195364
# PIN_1194_GC114074_SI-GA-E10_AACCTGAAGCTTAAGA-1  0.3968694 -0.8522779
# PIN_1194_GC114074_SI-GA-E10_AAGGAATTCCTATTGT-1 -2.0859216 -1.0033504
# PIN_1194_GC114074_SI-GA-E10_AAGTCGTTCTGAATGC-1 -1.4452787 -1.3332216
# PIN_1194_GC114074_SI-GA-E10_AATCGACTCTTTGATC-1 -3.2186957 -1.3799822


plot_data_srt=as.data.frame(SPATIAL[["umap"]]@cell.embeddings)
head(plot_data_srt)

#                       umap_1    umap_2
# PIN1561b_aaaahcno-1 6.664046 -6.025258
# PIN1561b_aaabdffa-1 1.406551  9.431500
# PIN1561b_aaachiil-1 2.009092  8.986633
# PIN1561b_aaacnkll-1 5.889764 -3.418880
# PIN1561b_aaadncjj-1 9.170426 -1.048979
# PIN1561b_aaaebebo-1 1.297517  9.165634

plot_data <- plot_data_srt
Col=plot_data$umap_1
Row=plot_data$umap_2

plan("multicore", workers = 15) 
options(future.globals.maxSize = 400000 * 1024^2) 

table(SC@meta.data$Cluster)

dot.srt <- setup.srt(srt_data = SPATIAL, srt_coords = plot_data, verbose=TRUE)
dot.ref <- setup.ref(ref_data = SC, ref_annotations = SC@meta.data$Cluster, max_genes = 500,ref_subcluster_size = 1,remove_mt = FALSE, verbose=TRUE) 
dot <- create.DOT(dot.srt, dot.ref)
dot <- run.DOT.highresolution(dot, iterations = 100,verbose=TRUE, ratios_weight = 0.75) # # Abundance weight; a larger value more closely matches the abundance of cell types in the spatial data to those in the reference data

dim(dot@weights)
#spatial cells x celltype
[1] 632242     12

dot@weights[0:3,0:5]
#                     B cells Endothelial cells Epithelial cells Fibroblasts
# PIN1561_aaaahcno-1 0.2726878      1.902432e-06     3.260303e-09  0.09741889
# PIN1561_aaabdffa-1 0.1917019      3.452548e-07     3.260303e-09  0.04538347
# PIN1561_aaachiil-1 0.2136111      1.291247e-06     3.260303e-09  0.06142429
#                            Glia
# PIN1561_aaaahcno-1 3.260303e-09
# PIN1561_aaabdffa-1 4.084539e-01
# PIN1561_aaachiil-1 2.104871e-01

saveRDS(dot,"dot_labels_ratio075.rds")
dot=readRDS("dot_labels_ratio075.rds")

library(tidyverse)
weights=as.data.frame(dot@weights)
head(weights)

plot_data[0:5,0:2]
#                      umap_1     umap_2
# PIN1561_aaaahcno-1 6.418548  5.0077626
# PIN1561_aaabdffa-1 2.808654 -9.9406785
# PIN1561_aaachiil-1 3.578207 -9.4499730
# PIN1561_aaacnkll-1 6.329224  2.9142004
# PIN1561_aaadncjj-1 8.858322  0.9711915

dim(plot_data)
dim(dot@weights)
dim(weights)
# same number of rows
# [1] 632242      2
# [1] 632242     12
# [1] 632242     12
head(rownames(plot_data))
head(rownames(weights))
# same and same order

df=merge(plot_data, weights, by="row.names", all=TRUE) 
head(df)
#            Row.names   umap_1     umap_2   B cells Endothelial cells
# 1 PIN1561_aaaahcno-1 6.418548  5.0077626 0.2726878      1.902432e-06
# 2 PIN1561_aaabdffa-1 2.808654 -9.9406785 0.1917019      3.452548e-07
# 3 PIN1561_aaachiil-1 3.578207 -9.4499730 0.2136111      1.291247e-06
# 4 PIN1561_aaacnkll-1 6.329224  2.9142004 0.1666628      4.437834e-07
# 5 PIN1561_aaadncjj-1 8.858322  0.9711915 0.2636895      1.134791e-02
# 6 PIN1561_aaaebebo-1 2.576019 -9.8102273 0.1228443      3.260303e-09
#   Epithelial cells Fibroblasts         Glia          ILC   Mast cells
# 1     3.260303e-09  0.09741889 3.260303e-09 1.212833e-02 1.760589e-02
# 2     3.260303e-09  0.04538347 4.084539e-01 4.655727e-02 3.260303e-09
# 3     3.260303e-09  0.06142429 2.104871e-01 2.398978e-07 5.007831e-02
# 4     4.710563e-09  0.12519702 3.260303e-09 2.738657e-03 1.447575e-02
# 5     7.042257e-03  0.04303707 4.694839e-03 3.012536e-02 3.286399e-02
# 6     3.260303e-09  0.05203469 4.057153e-01 1.032867e-01 3.260303e-09
#   Myeloid cells           NK Plasma cells Proliferating cells    T cells
# 1   0.074334902 7.042305e-03   0.22457015        1.241398e-07 0.29420971
# 2   0.167840379 5.219083e-08   0.04577465        1.105207e-07 0.09428795
# 3   0.112284823 8.115107e-08   0.12597809        5.536803e-08 0.22613459
# 4   0.109937405 6.651020e-03   0.33529069        1.760650e-02 0.22143975
# 5   0.104068861 1.877944e-02   0.22026622        4.303603e-03 0.25978091
# 6   0.001564948 2.269179e-02   0.14319249        3.260303e-09 0.14866980

head(rownames(df))
rownames(df)=df$Row.names
head(rownames(df))
df$Row.names=NULL
head(df)

df2=df %>% 
  rownames_to_column('id') %>%
  gather(celltype, score, 4:15) %>%  # columns number # ! start with rowname as 1
  group_by(id) %>% 
  filter(score == max(score)) %>% 
  arrange(id)
  
head(print(df2,width = Inf))

#   id                 umap_1 umap_2 celltype     score
#   <chr>               <dbl>  <dbl> <chr>        <dbl>
# 1 PIN1561_aaaahcno-1   6.42  5.01  T cells      0.294
# 2 PIN1561_aaabdffa-1   2.81 -9.94  Glia         0.408
# 3 PIN1561_aaachiil-1   3.58 -9.45  T cells      0.226
# 4 PIN1561_aaacnkll-1   6.33  2.91  Plasma cells 0.335
# 5 PIN1561_aaadncjj-1   8.86  0.971 B cells      0.264
# 6 PIN1561_aaaebebo-1   2.58 -9.81  Glia         0.406

min(df2$score)
[1] 0.1557121
max(df2$score)
[1] 0.9999997

dim(df2)
potential_numbers = df2 %>% group_by(celltype) %>%  summarise(n=n())  
potential_numbers 

df3= as.data.frame(df2)

# we have duplicate rows, see them
duplicates <- df3 |>
  group_by(id) |>
  filter(n() > 1) |>
  ungroup()

print(duplicates,n = Inf)

df4=df3 %>% distinct(id, .keep_all = TRUE)
dim(df4)
head(df4)

df3=df4

rownames(df3)=df3$id
df3$id=NULL
head(df3)

head(df3)
dim(df3)
write.table(df3, "spatial_cells_with_annotations.txt",col.names=NA, sep="\t")

# then only substract Mast cells - the other as "OTHERS"
df3$celltype[df3$celltype!="Mast cells"]= "other"
head(df3)
unique(df3$celltype)
write.table(df3, "spatial_cells_with_annotations_MAST_or_other.txt",col.names=NA, sep="\t")

plot_data=df3

# combine this dotr results with our Xenium.big

Col=plot_data$umap_1
Row=plot_data$umap_2

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "#E5E4E2", # platiunium - ligth grey
   "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown","black"
)

head(plot_data)
plot_data2 <- plot_data %>% arrange(desc(celltype))
head(plot_data2)
  
pdf("srt_coord_celltype.pdf", width=20, height=20)
ggplot(plot_data2, aes(x = Col, y = Row, color = celltype))+
  geom_point(size = 2)+
  theme_bw()+
  scale_colour_manual(values=c25)+
  guides(color = guide_legend(override.aes = list(size = 30)))+
  theme(panel.background = element_rect(fill = 'white'), 
        panel.grid = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.title=element_text(size=30),
        legend.text = element_text(size=30)
        )
dev.off()

# now, our spatial cells have embeddings + celltype annotation
# put them back on slides

Xenium.big=readRDS("/.../MERGED.rds")
names<-c("PIN1561",
"PIN1825b",
"PIN1913",
"PIN1916",
"PIN1861",
"PIN1635",
"PIN1914",
"PIN4063")
for (i in names){
  DefaultBoundary(Xenium.big[[i]]) <- "segmentation" #prefered over centroids
}


Xenium.big@meta.data$orig.ident=""
Xenium.big@meta.data$orig.ident=Xenium.big@meta.data$slide
head(Xenium.big@meta.data,2)
tail(Xenium.big@meta.data,2)

# add in metadata the celltype annotations you just calculated with DOT --> plot_data is the dataframe with the cell names and celltypes
head(plot_data)

# add to Xenium.big the plot_data info
Xenium_with_anno <- AddMetaData(object = Xenium.big, metadata = plot_data)
dim(Xenium_with_anno@meta.data)
head(Xenium_with_anno@meta.data,2)
table(Xenium_with_anno@meta.data$celltype)

DefaultAssay(Xenium_with_anno)="Xenium"
Xenium_with_anno

y=table(Xenium_with_anno@meta.data$celltype)
y
write.table(y, "cell_per_cluster_SPATIAL.txt",col.names=NA, sep="\t")

y=table(Xenium_with_anno@meta.data$celltype,Xenium_with_anno@meta.data$slide)
y
write.table(y, "cell_per_cluster_SPATIAL_per_slide.txt",col.names=NA, sep="\t")

saveRDS(Xenium_with_anno,"Xenium_with_anno_075.rds")
Xenium_with_anno[["PIN1561"]]
head(Xenium_with_anno@meta.data,2)

for (names in c("PIN1561","PIN1635","PIN1825b","PIN1861","PIN1913","PIN1914","PIN1916","PIN4063")){
pdf(paste0(names[1],"_ANNOTATIONS.pdf"), width=20, height=10)
print(ImageDimPlot(Xenium_with_anno, fov = names[1], group.by="celltype",border.size=NULL,border.color=NA, axes=TRUE)+
scale_fill_manual(name="Legend", values = c("Mast cells" = "#FF33FF","other"=NULL)))
dev.off()
}

# ========================= #
# NOW, WE CREATE SIGNATURES #
# ========================= #

library(DOTr)
library(ggplot2)

Xenium_with_anno=readRDS("Xenium_with_anno_075.rds")
Xenium.big=readRDS("/.../MERGED.rds")
dot=readRDS("dot_labels_ratio075.rds")
plot_data=read.table("spatial_cells_with_annotations.txt",header=TRUE,row.names = 1, sep="\t")

dot@weights[0:3,0:5]
#                      B cells Endothelial cells Epithelial cells Fibroblasts
# PIN1561_aaaahcno-1 0.2726878      1.902432e-06     3.260303e-09  0.09741889
# PIN1561_aaabdffa-1 0.1917019      3.452548e-07     3.260303e-09  0.04538347
# PIN1561_aaachiil-1 0.2136111      1.291247e-06     3.260303e-09  0.06142429
#                            Glia
# PIN1561_aaaahcno-1 3.260303e-09
# PIN1561_aaabdffa-1 4.084539e-01
# PIN1561_aaachiil-1 2.104871e-01

MC1_signature=c("CADPS","CCL20","CD83","CDK15","GATA2","GPR183","HDC","IL1B","IL1RAPL1","MLPH","PTTG1","THBS1")
MC2_signature=c("CADPS","CCL20","CD83","CDK15","GATA2","GPR183","HDC","IL1B","IL1RAPL1","MLPH","PTTG1","IER5")
MC3_signature=c("AREG","BATF","CALB2","CCL4","CMA1","CTSG","RPS4Y1","ID2","SOCS3")
MC4_signature=c("AREG","BATF","CALB2","CMA1","CTSG","RPS4Y1","CST7")
MC5_signature=c("APOE")

Xenium.big=AddModuleScore(Xenium.big,features=list(c(MC1_signature)),ctrl = 5,name = "MC1_signature")
names(Xenium.big@meta.data)[names(Xenium.big@meta.data) == "MC1_signature1"] <- "MC1_signature"
Xenium.big=AddModuleScore(Xenium.big,features=list(c(MC2_signature)),ctrl = 5,name = "MC2_signature")
names(Xenium.big@meta.data)[names(Xenium.big@meta.data) == "MC2_signature1"] <- "MC2_signature"
Xenium.big=AddModuleScore(Xenium.big,features=list(c(MC3_signature)),ctrl = 5,name = "MC3_signature")
names(Xenium.big@meta.data)[names(Xenium.big@meta.data) == "MC3_signature1"] <- "MC3_signature"
Xenium.big=AddModuleScore(Xenium.big,features=list(c(MC4_signature)),ctrl = 5,name = "MC4_signature")
names(Xenium.big@meta.data)[names(Xenium.big@meta.data) == "MC4_signature1"] <- "MC4_signature"
Xenium.big=AddModuleScore(Xenium.big,features=list(c(MC5_signature)),ctrl = 5,name = "MC5_signature")
names(Xenium.big@meta.data)[names(Xenium.big@meta.data) == "MC5_signature1"] <- "MC5_signature"
head(Xenium.big@meta.data)

Xenium_with_anno <- AddMetaData(object = Xenium.big, metadata = plot_data)

DefaultAssay(Xenium_with_anno)="Xenium"
Xenium_with_anno

Idents(Xenium_with_anno)="celltype"
Xenium_mast=subset(Xenium_with_anno,idents="Mast cells")
Xenium_mast

df=Xenium_mast@meta.data[,16:22]
head(df,5)

#                    MC1_signature MC2_signature MC3_signature MC4_signature
# PIN1561_aacameaa-1   -0.08223780   -0.07048954    0.07526598   -0.02038668
# PIN1561_aacnpdkp-1   -0.04996837   -0.09598238    0.02054120    0.05824766
# PIN1561_abafpnlp-1   -0.09696140   -0.08212295    0.25584279    0.29662860
# PIN1561_abahfbif-1   -0.10285467   -0.10663677    0.19954237    0.24436427
# PIN1561_abchofbl-1    0.01775638    0.05300115   -0.04726004   -0.02038668
#                    MC5_signature    umap_1    umap_2
# PIN1561_aacameaa-1    -0.1386294 10.047413 -4.336563
# PIN1561_aacnpdkp-1    -0.1386294 10.335008 -4.275501
# PIN1561_abafpnlp-1    -0.1386294  4.148210 -6.701655
# PIN1561_abahfbif-1     0.9599829  4.161855 -6.734832
# PIN1561_abchofbl-1    -0.3583519 10.232499 -4.276386

df2=df %>% 
  rownames_to_column('id') %>%
  gather(mast_subpop, score, 2:6) %>%  # columns number # ! start with rowname as 1
  group_by(id) %>% 
  filter(score == max(score)) %>% 
  arrange(id)
  
# we have duplicate rows, see them
duplicates <- df2 |>
  group_by(id) |>
  filter(n() > 1) |>
  ungroup()

print(duplicates,n = Inf)

df4=df2 %>% distinct(id, .keep_all = TRUE)
dim(df4)
head(df4)
df2=df4


xx=dplyr::filter(df2, score >0)
dim(xx)

potential_numbers = xx %>% group_by(mast_subpop) %>%  summarise(n=n())  
potential_numbers 
write.table(xx, "cells_with_anno_C1_5_SMALL_SIGNATURE.txt",col.names=NA, sep="\t")
plot_data=xx

# first, I need plot_data to be the same length as the seurat object
# make a new dataframe
head(Xenium.big@meta.data,2)
df_to_merge=Xenium.big@meta.data[,1:3]
head(df_to_merge)
df_to_merge2 <- tibble::rownames_to_column(df_to_merge, "id") 
head(df_to_merge2)


# merge plot_data with df_to_merge
df_to_merge3=as.data.frame(plot_data)
head(df_to_merge3)


library(dplyr)
DF=left_join(df_to_merge2, df_to_merge3 , by = "id")
DF$mast_subpop[is.na(DF$mast_subpop)] <- "other"

head(DF,30)
yy=table(DF$mast_subpop)

rownames(DF)=DF$id
DF$id=NULL
head(DF)

DF2=DF[,c(4:7)]
plot_data=DF2
head(plot_data)
# add to Xenium.big the plot_data info
Xenium_with_anno <- AddMetaData(object = Xenium.big, metadata = plot_data)
dim(Xenium_with_anno@meta.data)
# [1] 632242     24
head(Xenium_with_anno@meta.data,2)
table(Xenium_with_anno@meta.data$mast_subpop)


Col=plot_data$umap_1
Row=plot_data$umap_2

c25 <- c(
  "dodgerblue2",  # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
   # platiunium - ligth grey
   "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown","black"
)

head(plot_data)
plot_data2 <- plot_data %>% arrange(desc(mast_subpop))
head(plot_data2)
  
pdf("srt_coord_celltype_MCsubsets.pdf", width=20, height=20)
ggplot(plot_data2, aes(x = Col, y = Row, color = mast_subpop))+
  geom_point(size = 2)+
  theme_bw()+
  scale_colour_manual(values=c25)+
  guides(color = guide_legend(override.aes = list(size = 30)))+
  theme(panel.background = element_rect(fill = 'white'), 
        panel.grid = element_blank(),
        axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        legend.title=element_text(size=30),
        legend.text = element_text(size=30)
        )
dev.off()

DefaultAssay(Xenium_with_anno)="Xenium"
Xenium_with_anno
y=table(Xenium_with_anno@meta.data$mast_subpop)
y
write.table(y, "cell_per_cluster_SPATIAL.txt",col.names=NA, sep="\t")
y=table(Xenium_with_anno@meta.data$mast_subpop,Xenium_with_anno@meta.data$slide)
y
write.table(y, "cell_per_cluster_SPATIAL_per_slide.txt",col.names=NA, sep="\t")
saveRDS(Xenium_with_anno,"Xenium_with_anno_SUBPOP_075_Sign.rds")
head(Xenium_with_anno@meta.data,2)

for (names in c("PIN1561","PIN1635","PIN1825b","PIN1861","PIN1913","PIN1914","PIN1916","PIN4063")){
pdf(paste0(names[1],"_ANNO_SUBPOP_SMALL_Signatures.pdf"), width=20, height=10)
print(ImageDimPlot(Xenium_with_anno, fov = names[1], group.by="mast_subpop",border.size=NULL,border.color=NA, axes=TRUE)+
scale_fill_manual(name="Legend", values = c("MC1_signature" = "#00CC66","MC2_signature" = "#CCFF33","MC3_signature" = "#FF0000","MC4_signature"="#FF33FF","MC5_signature"="#0000FF","other"=NULL)))
dev.off()
}
