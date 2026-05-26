# 1) Gprofiler2

suppressPackageStartupMessages({library(Seurat)
 suppressWarnings(library(ggplot2))
 suppressWarnings(library(stringr))
 suppressWarnings(library(ggrepel))
 suppressWarnings(library(cowplot))
 suppressWarnings(library(SingleCellExperiment))
 suppressWarnings(library(dplyr))
 suppressWarnings(library(patchwork))
 suppressWarnings(library("readxl"))
 suppressWarnings(library(gprofiler2))   # v0.2.3
 suppressWarnings(library(ggpubr))                 
 })

setwd("/path") 
Seurat_Final=readRDS("OBJECT.rds")
genes=read.table("wilcox_markers_MC1toMC5.txt", sep="\t", header=TRUE)

setwd("/path/GProfiler") 
sub=Seurat_Final

genes2 <- genes[genes$p_val_adj<0.05,]
genes3 <- genes2[genes2$avg_log2FC >0,]
genes3 <- genes3[order(genes3$avg_log2FC, decreasing = T),]
head(genes3)
clusters_sub<-levels(sub$Cluster)

SELECTED=c("NFKB1-NFKB2-REL-RELA-RELB complex","programmed cell death","response to topologically incorrect protein","response to cytokine","antigen processing and presentation","Unfolded Protein Response (UPR)","positive regulation of nitrogen compound metabolic process","positive regulation of RNA metabolic process","negative regulation of transcription by RNA polymerase II","Activation of the AP-1 family of transcription factors","response to stress","granulocyte chemotaxis","canonical NF-kappaB signal transduction","regulation of cell-cell adhesion","regulation of cysteine-type endopeptidase activity","apoptotic process","response to tumor necrosis factor","monocyte chemotaxis","angiotensin maturation","positive regulation of calcium ion import","eosinophil chemotaxis","blood vessel morphogenesis","chylomicron remnant clearance","very-low-density lipoprotein particle clearance","triglyceride-rich lipoprotein particle clearance","positive regulation of heparan sulfate proteoglycan binding","Cholesterol metabolism")

go.list<-list()
SOURCE="ALL"
for(j in 1:length(clusters_sub)){
    tryCatch({
   print(j)
   xxx=clusters_sub[j]
   print(xxx)    
     Cluster_go<-genes3[which(genes3$cluster %in% clusters_sub[j] ),]$gene 
     gostCluster_go<-gost(query=Cluster_go, organism = "hsapiens", sources=c("REAC","GO:BP","CORUM","KEGG"), evcodes = T,ordered_query=T)     
     df_go<- as.data.frame(gostCluster_go$result)
     df_go<- apply(df_go,2,as.character)
     go_Cluster_results_up<-as.data.frame(gostCluster_go$result)
     go_Cluster_results_up$p_value_log<- -log10(go_Cluster_results_up$p_value) 
     go_Cluster_results_up = filter(go_Cluster_results_up, term_name %in%  SELECTED)
  
     df2 <- go_Cluster_results_up[order(go_Cluster_results_up$intersection_size,decreasing=FALSE),]
     head(df2)
     tail(df2)
     go_Cluster_results_up=df2
     up_plot<-ggbarplot(go_Cluster_results_up, x = "term_name", y = "intersection_size",fill=c("steelblue2"),label = F, label.pos = "out",orientation = "horiz", title=paste0(clusters_sub[j], " selected pathways", sep=""))
     go.list[[j]]<-up_plot
}   ,error=function(e){})}  
  
plot_cluster_go<-patchwork::wrap_plots(plots = go.list, nrow=3)
ggsave(plot=plot_cluster_go, filename=paste0("MC_HV_Rectum_Colon_gProfiler_selected_intersection_",SOURCE,".png") ,height=20, width=15*(length(clusters_sub)/3), units="in", dpi=320)    
ggsave(plot=plot_cluster_go, filename=paste0("MC_HV_Rectum_Colon_gProfiler_selected_intersection_",SOURCE,".pdf") ,height=20, width=15*(length(clusters_sub)/3), units="in", dpi=320)    

# 2) UCell

library(UCell)

gene.sets <- list(MC1=c("HSPA5","PRG2","FAM107B","LMNA","HSP90B1","BCL2A1","HERPUD1","RRAD","SDF2L1","THBS1","NFKB1","TUBA1A","PIM3","IL13","RGCC","REL","NFE2L2","IGHA1","CDC42EP3","MAPK6","GLUL","SRSF7","HES4","PDIA3","CALR","HLA-B","MYADM","KDM6B","SKIL","TNFAIP8","EFHD2","RELB","DUSP6","DNAJB11","XBP1","HSPH1","AHR","TSC22D1","HLA-A","SDCBP","KMT2E","NINJ1","VIM","JUND","HDC","GALNT6","YBX3","NFKB2","HECW2-AS1","TNFRSF9"),
MC2=c("JUN","FOS","DUSP1","KLF2","RGS1","HSPA1B","EGR1","KLF4","DNAJB1","FOSB","PPP1R15A","HSPA1A","JUNB","IER2","TSC22D3","ZFP36","SEMA4A","INTS6","AC020916.1","BTG2","HES1","SLC38A2","NEU1","MTRNR2L12","UBE2S","SERTAD1","KLF6","CD69","BTG1","PPP1R10","IRF1","KCNQ1OT1","IGHA1","CSF1","SOX4","DNAJA1","HEXIM1","ATF3","NR4A2","NR4A1","ZFP36L2","DDX3X","DUSP2","MTRNR2L8","SAT1","TUBB4B","TXNIP","EGR3","CITED2","BRD2"),
MC3=c("CCL4","CCL4L2","MT2A","AREG","ID2","LGALS1","BIRC3","BATF","ARID5B","VEGFA","PLAUR","S100A10","CCL2","NTM","S100A4","PTGS2","CD52","ATP1B3","RHOH","IFI16","IL2RA","SRGN","CPEB4","TMOD1","LGALS3","S100A6","AL445524.1","CALB2","SOCS3","FTH1","AP1S3","G0S2","RNF145","C5AR1","CYTOR","GABPB1","IFITM3","SPATA13","EMP3","PAG1","NAMPT","EGLN1","CREM","REL","CCNH","L1CAM","KLF7","HINT1","BHLHE40","C21orf91"),
MC4=c("CMA1","HPGD","CCL2","LGALS3","AC253572.2","CKLF","GSTM3","GLRX","NTM","CTSG","AREG","HINT1","MT2A","SNU13","GABARAP","NDUFA4","PEBP1","ATP5MC2","EVI2A","CHCHD10","AL078590.2","ABRACL","TESC","ASAH1","PLIN2","AC068587.4","OST4","IER2","DAD1","TYROBP","COX7A1","FOS","ZFP36L2","CST7","SRGN","CDCA7L","TPT1","UQCRQ","AL355881.1","NDUFC1","FCER1G","RGS10","FDXR","CHCHD6","CSTB","NPC2","NDUFA3","ATP6V1G1","DNAJC19","TOMM7"),
MC5=c("APOE","CD52","LGALS1","APOC1","PMEPA1","BCL2A1","HAS2","LTC4S","RNASE1","S100A4","FNBP1","MALT1","NCOA4","NSMCE1","PLAUR","CST3","S100A11","ARMH1","PPP1R14B","S100A10","TMSB4X","THAP2","LAPTM5","ATP5F1E","LIMS1","TNFRSF18","GAPDH","HES4","IFI30","SH3BGRL3","CORO1A","LAT","CLU","TYROBP","ARL6IP5","TLN1","SQLE","PFN1","TSPO","RGS10","GMFG","CAVIN1","S100A6","APOC2","TNFRSF4","RNASEK","CD99","RHOG","DRAP1","TMSB10")
)

UCell <- AddModuleScore_UCell(Seurat_Final, features = gene.sets)
signature.names <- paste0(names(gene.sets), "_uCell")
head(UCell@meta.data,2)

DF2= UCell@meta.data
head(DF2,2)

target <- c("MC1","MC3","MC4","MC5","MC2")
library(dplyr)
DF=DF2 %>% arrange(factor(Cluster, levels = target))
head(DF,10)
tail(DF,20)

X_axis=DF$MC1_UCell
Y_axis=DF$MC3_UCell

title_X="MC1_UCell"
title_Y="MC3_UCell"

pdf(paste0("TOP50_scatter_with_",title_Y,"_vs_",title_X,".pdf"), width=8, height=8)
ggplot(DF, aes(x = X_axis, y = Y_axis,colour = factor(Cluster))) +
geom_point() + 
xlab(paste0(title_X," (%)")) + ylab(paste0(title_Y," (%)")) +
scale_color_manual(values=c("#F8766D","#F5B700","#00A9FF","#00BE67","#ED68ED")) 
dev.off()

# 3) OTHER PLOTS

# e.g. feature plot

for (gene in c("KIT","PTGS1","HPGDS")){
    tryCatch({
pdf(paste0("Specific_gene_",gene,"_Feature.pdf"),width = 8, height = 8)
print(SCpubr::do_FeaturePlot(Seurat_Final, features = gene, reduction = "umap", pt.size = 2, order=TRUE,raster = FALSE)+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")) ))#,limits = c(0,3)))
dev.off()
},error=function(e){})
}

# e.g. Dotplot

library("ggplot2")
library(SCpubr)

GROUP="xxx"
genes=list(
    "Cytokine receptors"=(c("IL1RL1","IL2RA","IL2RG","IL4R","IL5RA","IL6R","IL6ST","IL10RB","IL18R1","IL27RA","IFNGR1","IFNGR2","IFNAR1","IFNAR2","TGFBR1","TGFBR2")),
    "Ag presentation"=(c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-A","HLA-B","HLA-C","HLA-E")),
    "Co-stimulation"=(c("CD80","CD86","CD70","CD40","CD48","CD58","ICAM1","CD83","CD274")),
    "Cytoskeleton"=(c("TUBA1A","TUBA1B","TUBB4B","TUBB6","ACTG1","CDC45EP3","VASP","SAMSN1","EFHD2","LMNA","TAGLN2"))    
)

Idents(Seurat_Final)="Cluster"
my_levels <- c("MC5","MC4","MC3","MC2","MC1")

# Re-level object@ident
Seurat_Final@active.ident <- factor(x = Seurat_Final@active.ident, levels = my_levels)

p1 <- DotPlot(Seurat_Final, features = (genes), dot.scale=10, cols="RdBu",scale=TRUE) + theme(text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
p2=p1+theme(strip.background = element_rect(color="black", linetype="solid"), axis.title= element_blank()) 
p3=p2+ scale_color_distiller(direction=-1, palette="RdBu")

pdf(paste0("MC_HV_Rectum_Colon_",GROUP,"_Dotplot_per_cluster.pdf"), height=3.5, width=12+5)
p3
dev.off()

# e.g. Violin plot

# subset MC5 only
Idents(Seurat_Final)="Cluster"
sub=subset(Seurat_Final, idents="MC5")
sub 
Idents(sub)="Layer"
genes=c("APOE","CD52","LGALS1","APOC1","GPR35")

for (i in (1:length(genes))) {                                                                                               
gene=genes[i]     
tryCatch({
pdf(paste0("MC_HV_Rectum_Colon_",gene,"_Vln_per_LAYER_inMC5only.pdf"),width = 6, height = 6)
print(VlnPlot(sub, features = gene, pt.size = 0) + theme(axis.text.x = element_text(angle = 90, size=10))
)
dev.off()
},error=function(e){})}

