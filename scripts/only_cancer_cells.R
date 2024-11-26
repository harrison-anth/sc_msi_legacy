#generating pseudo-bulk
#functions and libs

library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

set.seed(seed = 152727)


process_data <- function(sample_name,filter_norm){
  print(paste("Begin processing sample",sample_name,sep=" "))
  print(paste("Read in file",sample_name,sep=" "))
  temp <- Read10X_h5(filename = paste0('/data4/hanthony/single_msi/cell_ranger_output/',sample_name,
'/outs/filtered_feature_bc_matrix.h5'))
  print(paste("Converting to Seurat object, filtering, and normalizing"))
  temp <- CreateSeuratObject(counts=temp,
                             project="MSI",
                             min.cells=3,
                             min.features=200)
  #add percent mitochondrial
  temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern = "^MT-")
  #add atomic_calls
  atomic_cells <- readRDS(paste0('../atomic/',sample_name,'.rds'))
  temp <- AddMetaData(temp, metadata=atomic_cells$pan_cancer_cluster,col.name='tumor')
  temp <- subset(temp,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

  if(filter_norm == TRUE){

  temp <- subset(temp,subset=tumor == "Cancer")
	}

  temp <- NormalizeData(temp,normalization.method = "LogNormalize",scale.factor = 10000)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
  
  #plot most variable genes
  #top10 <- head(VariableFeatures(temp),10)
  #LabelPoints(plot=(VariableFeaturePlot(temp),points=top10,repel=TRUE))
  #scale data
  all.genes <- rownames(temp)
  #might still need to "regress out" mitchonodrial and other sources of variation
  temp <- ScaleData(temp,feature=all.genes)
  temp <- RunPCA(temp, features=VariableFeatures(object=temp),npcs=5)
  
  
  print(paste("Identifying clusters"))
  
  temp <- FindNeighbors(temp,dims=1:5)
  temp <- FindClusters(temp,resolution = 0.5)
  temp <- RunUMAP(temp, dims = 1:5)
  #find marker genes that are differentially expressed in clusters
   
  
  
  
  
  
  return(temp)
  
  
}

add_sensor_rna <- function(s_obj,sample_name,threshold){

sensor_rna_results <- fread(paste0('../sensor_rna_results/',sample_name,'.txt'))

sensor_rna_results$msi_status[sensor_rna_results$`probability_of_MSI-H` >= threshold ] <- 'MSI-H'
sensor_rna_results$msi_status[sensor_rna_results$`probability_of_MSI-H` < threshold ] <- 'MSS'



temp <- t(sensor_rna_results)
colnames(temp) <- temp[1,]
temp <- temp[2:3,]
s_obj <- AddMetaData(s_obj,temp[1,],col.name='sensor_rna_status')
s_obj <- AddMetaData(s_obj,temp[2,],col.name='sensor_rna_prob')

saveRDS(s_obj, paste0('../annotated_h5/',sample_name,'_cancer_only.rds'))


return(s_obj)

}



new_sample <- process_data(sample_name,filter_norm=TRUE)

new_sample2 <- add_sensor_rna(new_sample,sample_name,threshold=.75)


