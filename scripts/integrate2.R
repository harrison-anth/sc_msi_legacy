#Integrate sample from common patient ID's
#Much of the integration code from user alg697's comment on seurat github issue
#https://github.com/satijalab/seurat/issues/8938

#functions and libs
library(scATOMIC)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(R.utils)
library(scATOMIC)
library(tidyverse)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
gsm <- as.character(argus[2])





set.seed(seed = 152727)

if(gsm == 'Y'){key <- fread('../manifests/final_gsm_key.tsv',header=TRUE)
} else if(gsm == 'N'){key <- fread('../manifests/final_key.tsv',header=TRUE)}


integrate_data <- function(sample_name){

indv_key <- filter(key, patient_id == sample_name)
  print(paste("Begin processing sample",sample_name,sep=" "))
  print(paste("Integrating", nrow(indv_key),"samples", sep=" "))
if(nrow(indv_key)<2){
stop('Only 1 sample no need to integrate')
}
all_objs <- list()
 
for(i in 1:nrow(indv_key)) {

if(gsm == 'Y'){
file_name <- indv_key$sample_id[i]
} else if (gsm == 'N'){file_name <- indv_key$filename[i]}



assign(x = paste0('s_obj'),value = readRDS(paste0('../annotated_h5/',file_name,'.rds')))

atomic_cells <- readRDS(paste0('../atomic/',file_name,'.rds'))

s_obj$scATOMIC_pred <- atomic_cells$scATOMIC_pred

s_obj$classification_confidence <- atomic_cells$classification_confidence


if('pan_cancer_cluster' %in% colnames(atomic_cells@meta.data)){
s_obj$pan_cancer_cluster <- atomic_cells$pan_cancer_cluster
}


s_obj$sensor_rna_prob <- as.numeric(s_obj$sensor_rna_prob)

assign(x = paste0(sample_name,'_', i),value = s_obj) %>% append(all_objs) -> all_objs
}

merged <- merge(all_objs[[1]],
	y=c(all_objs[2:length(all_objs)]),
	 project = sample_name)
pre <- merged %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA(.) %>%
FindNeighbors(.,dims=1:10, reduction = "pca") %>% FindClusters(., resolution = c(0.5)) %>%
RunUMAP(.,dims=1:10, reduction = "pca")

#sample "CRC2817" could only work with k.weight=50 (as opposed to 100 default) due to low number of cells in some sample(s)

integrated <- IntegrateLayers(object = pre, method = CCAIntegration, orig.reduction = "pca", 
new.reduction="integrated.cca", verbose = TRUE,k.weight=50)

integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- integrated %>% FindNeighbors(., reduction = "integrated.cca", dims = 1:10) %>% 
FindClusters(.,res=c(0.5), cluster.name = "integrated.clusters") %>% 
RunUMAP(., dims=1:10, reduction= "integrated.cca", reduction.name = "integrated.umap")

saveRDS(integrated, paste0('../integrated_samples/',sample_name,'.rds'))

#now just the cancer
#recluster
cancer_int <- subset(integrated,subset=pan_cancer_cluster == "Cancer")
DefaultAssay(cancer_int) <- "RNA"
cancer_fin <- cancer_int %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
#ScaleData(., vars.to.regress = c("percent.mt")) %>%
ScaleData(.) %>%
RunPCA(.) %>% RunUMAP(., dims = 1:30) %>% FindNeighbors(., dims = 1:30) %>%
FindClusters(., resolution = 0.8)

saveRDS(cancer_fin, paste0('../integrated_samples/',sample_name,'_cancer.rds'))

}

integrate_data(sample_name)

