library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(ggplot2)
library(R.utils)
library(foreach)
library(parallel)
library(copykat)
library(cutoff.scATOMIC)
library(scATOMIC)
library(devtools)

set.seed(seed = 152727)

n_cores=detectCores()-1

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

reticulate::use_python("/data3/hanthony/.conda_envs/atomic/bin/python3",required=T)
sys <- import("sys")

rds1 <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))


sparse_matrix <- rds1@assays$RNA$counts


cell_predictions <- run_scATOMIC(sparse_matrix,mc.cores = n_cores)
results_temp <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = F, modify_results = T, mc.cores = n_cores, raw_counts = sparse_matrix, min_prop = 0.5 )

if( ! ('pan_cancer_cluster' %in% colnames(results_temp))){results_temp$pan_cancer_cluster <- 'Normal'}




rds2<- AddMetaData(rds1, results_temp)

saveRDS(rds2,paste0('../atomic/',sample_name,'.rds'))

#create cancer only barcodes

barcodes <- data.frame(barcode=names(rds2$pan_cancer_cluster),call=rds2$pan_cancer_cluster)

fwrite(barcodes,sep='\t',paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_cancer_barcodes.tsv'),row.names=FALSE,col.names=FALSE)


