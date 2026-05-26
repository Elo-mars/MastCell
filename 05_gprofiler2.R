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
