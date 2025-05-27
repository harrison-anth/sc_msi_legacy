library(Seurat)
library(PreMSIm)
library(janitor)
library(data.table)
library(tidyverse)
library(R.utils)


argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

path = system.file("extdata", "example.txt", package = "PreMSIm", mustWork = TRUE)
input_data <- data_pre(path, type = "ID")

#create count mat from s_obj
s_obj <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))
ct_mat <- t(as.matrix(s_obj@assays$RNA$counts))

#find premsim results for genes found in both training data and gene count matrix
msi_results <- msi_pre(ct_mat[,c(colnames(input_data)[colnames(input_data) %in% colnames(ct_mat)])])

prob <- attributes(msi_results$MSI_status)

s_obj$premsim_prob <- prob$prob

attr(msi_results$MSI_status, "prob") <- NULL
s_obj$premsim_status <- msi_results$MSI_status

saveRDS(s_obj, paste0('../premsim/',sample_name,'_premsim.rds'))


#for sensor_rna

ct_mat <- as.data.frame(ct_mat) %>% rownames_to_column('SampleID')

#create custom training data columns for this sample

pre_model <- fread ('/data4/hanthony/tcga_msi_tools/baselines/sensor_rna_files/model/demo/train.csv')

common_names <- intersect(colnames(ct_mat),colnames(pre_model))

new_model <- pre_model %>% select(all_of(c('SampleID','msi',common_names)))

fwrite(new_model,paste0('../temp/',sample_name,'_training_model.csv'),sep=',')

new_data <- ct_mat %>% select(all_of(c('SampleID',common_names)))

new_data2 <- new_data[ rowSums(new_data[,2:ncol(new_data)]) >0, ]

fwrite(new_data2,paste0('../temp/',sample_name,'.csv'),sep=',')

