library(Seurat)
library(janitor)
library(data.table)
library(tidyverse)
library(R.utils)


argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

#create count mat from s_obj
s_obj <- readRDS(paste0('../integrated_samples/',sample_name,'.rds'))
s_obj$All <- "1"
s_obj <- AggregateExpression(object = s_obj, group.by='All')
colnames(s_obj$RNA) <- sample_name
ct_mat <- t(as.matrix(s_obj$RNA))

#coerce rownames to column for sensorrna
ct_mat <- as.data.frame(ct_mat) %>% rownames_to_column('SampleID')

#create custom training data columns for this sample

pre_model <- fread ('../baselines/sensor_rna.csv')

common_names <- intersect(colnames(ct_mat),colnames(pre_model))

new_model <- pre_model %>% select(all_of(c('SampleID','msi',common_names)))

fwrite(new_model,paste0('../temp/',sample_name,'_glob_training_model.csv'),sep=',')

new_data <- ct_mat %>% select(all_of(c('SampleID',common_names)))

new_data2 <- new_data[ rowSums(new_data[,2:ncol(new_data)]) >0, ]

fwrite(new_data2,paste0('../temp/',sample_name,'_glob.csv'),sep=',')

