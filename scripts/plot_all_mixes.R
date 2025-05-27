library(Seurat)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(cowplot)
library(Nebulosa)
library(ggpubr)
library(R.utils)
library(RColorBrewer)
library(gridExtra)
set.seed(seed = 152727)


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,1))


mixing_table <- fread('../manifests/mixing_table_manifest.tsv')
mixing_table$file_prefix <- NULL

mix_names <- list.files(path="../annotated_h5/",pattern="*mix*",full.names=FALSE,ignore.case=TRUE)
mix_names <- gsub(mix_names,pattern=".rds",replacement="")

get_sens <- function(mix_name){
s_obj <- readRDS(paste0("../annotated_h5/",mix_name,".rds"))
print(FeaturePlot(s_obj,features=c('sensor_rna_prob'),cols=c('blue','red')) +
 ggtitle(paste(mix_name,' MSI score')))+sc
}
get_canc <- function(mix_name){
s_obj <- readRDS(paste0("../annotated_h5/",mix_name,".rds"))
print(DimPlot(s_obj,group.by='pan_cancer_cluster',cols=c('red','blue','grey'))+ggtitle(paste(mix_name,' ATOMIC')))
}
get_origin <- function(mix_name){
s_obj <- readRDS(paste0("../annotated_h5/",mix_name,".rds"))
s_obj2 <- readRDS(paste0("../temp/",mix_name,"_ground_truth.rds"))

s_obj$ground_truth <- s_obj2$ground_truth

print(DimPlot(s_obj,group.by='ground_truth')+ggtitle(paste(mix_name,' Origin')))
}

print("plotting MSI scores")
glist1 <- lapply(mix_names,get_sens)
print("plotting cancer status")
glist2 <- lapply(mix_names,get_canc)
print("plotting ground truth")
glist3 <- lapply(mix_names,get_origin)

pdf("../images/all_mixes.pdf",onefile=TRUE,height=12,width=26)

grid.table(mixing_table)

for(i in 1:length(glist1)){
print(marrangeGrob(grobs=c(glist1[i],glist2[i],glist3[i]),nrow=1,ncol=3))
}
#marrangeGrob(grobs=c(glist1,glist2,glist3),ncol=2,nrow=1)

dev.off()

#ggsave("../images/all_mixes.pdf",grid.arrange(grobs=c(glist1,glist2,glist3),ncol=3))
