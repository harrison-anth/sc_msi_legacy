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
gsm_key <- fread('../manifests/final_gsm_key.tsv')



set.seed(seed = 152727)


process_data <- function(sample_name,filter_norm,gsm){
  print(paste("Begin processing sample",sample_name,sep=" "))
  print(paste("Read in file",sample_name,sep=" "))


if(gsm==TRUE){
indv_key <- filter(gsm_key, sample_id == sample_name)
temp <- ReadMtx(
  mtx=paste0('../gsm_samps/',indv_key$file_prefix,'_matrix.mtx.gz'),
  cells=paste0('../gsm_samps/',indv_key$file_prefix,'_barcodes.tsv.gz'),
  features=paste0('../gsm_samps/',indv_key$file_prefix,'_features.tsv.gz')
)

} else if(gsm == FALSE){
  temp <- Read10X_h5(filename = paste0('/data4/hanthony/single_msi/cell_ranger_output/',sample_name,
'/outs/filtered_feature_bc_matrix.h5'))
 # print(paste("Converting to Seurat object, filtering, and normalizing"))
 } else{print('SAMPLE NAME ERROR')
stop()}
  temp <- CreateSeuratObject(counts=temp,
                             project="MSI",
                             min.cells=3,
                             min.features=100)




  #add percent mitochondrial
  temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern = "^MT-")
  #filter by minimum/maximum number of genes and percent mitochondrial
  #show pre-filter graphs
  #VlnPlot(temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #show feature scatterplot between any two features
  #FeatureScatter(temp,feature1="nCount_RNA",feature2="percent.mt")
  #print(paste('Filtering normal cells'))
  #copykat <- fread(paste0('copy_kat_results/',sample_name,'_copykat_prediction.txt'))
  #temp[['tumor']] <- NA
  #for (i in 1:length(names(temp$tumor))){
  #  temp2 <- filter(copykat, cell.names == names(temp$tumor[i]))
  #  temp$tumor[i] <- temp2$copykat.pred
  #  
  #}
  if(filter_norm==TRUE)
  {
    temp <- subset(temp,subset=tumor == "aneuploid")
  }
  temp <- subset(temp,subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 35)
  temp <- NormalizeData(temp,normalization.method = "LogNormalize",scale.factor = 10000)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
  
  #plot most variable genes
  #top10 <- head(VariableFeatures(temp),10)
  #LabelPoints(plot=(VariableFeaturePlot(temp),points=top10,repel=TRUE))
  #scale data
  all.genes <- rownames(temp)
  #might still need to "regress out" mitchonodrial and other sources of variation
  temp <- ScaleData(temp,feature=all.genes)
  temp <- RunPCA(temp, features=VariableFeatures(object=temp))
  
  
  print(paste("Identifying clusters"))
  
  temp <- FindNeighbors(temp,dims=1:15)
  temp <- FindClusters(temp,resolution = 0.5)
  temp <- RunUMAP(temp, dims = 1:15)
  #find marker genes that are differentially expressed in clusters
  
  
  
  
  
saveRDS(temp, paste0('../filtered_h5/',sample_name,'.rds'))
  
  return(temp)
  
  
}

#write out processed h5 as a new file for local review

gen_barcodes <- function(sample){

#supply an object returned by the process_data() function

barcodes <- names(sample$nCount_RNA)
clusters <- sample$seurat_clusters
num_clust <- length(levels(clusters)) -1
tempy <- data.frame(barcodes=barcodes,cluster=as.character(clusters))
fwrite(tempy,file=paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'),
sep='\t',col.names=FALSE,row.names=FALSE)

  for(z in 0:num_clust){
    temp <- as.matrix(t(clusters))
    barcodes2 <- names(temp[,temp[] == z])
    fwrite(data.table(barcodes2),file = paste0('../pseudobulk_barcodes/',sample_name,
	'/',sample_name,'_cluster_',z,'.tsv'),sep='\t',col.names = FALSE)
    }
}


<<<<<<< HEAD
#determine if sample is gsm or sra
if(sample_name %in% gsm_key$sample_id){
gsm=TRUE} else{gsm=FALSE}

tryme <- process_data(sample_name,FALSE,gsm)
=======
tryme <- process_data(sample_name,FALSE)
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef

gen_barcodes(tryme)
