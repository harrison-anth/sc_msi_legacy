#annotate filtered h5 files
#load libraries
library(numbat)
library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)
library(data.table)
library(maftools)
library(foreach)
library(parallel)

#add command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

#set seed
set.seed(seed = 152727)

#functions to annotate h5 for each MSI tool and ploidy results
add_sensor2 <- function(s_obj,sample_name){
  s_obj[['sensor2']] <- NA
    cluster_info <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'))
    cluster_info <- data.frame(clusters = sort(unique(cluster_info$V2)))
        for(i in 1:nrow(cluster_info)){
   num <- cluster_info$clusters[i]
      temp15 <- fread(paste0('../sensor2_results/',sample_name,'_cluster_',num))
      cluster_info$msi[i] <- temp15$`%`
      cluster_info$tot_sites[i] <- temp15$Total_Number_of_Sites
	if(cluster_info$tot_sites[i] < 10 ){cluster_info$msi[i] <- NA}
}
for(i in 1:nrow(s_obj[['sensor2']])){
  temp_name <- names(s_obj$sensor2[i])
  temp_cluster <- as.numeric(as.character(s_obj$seurat_clusters[temp_name == names(s_obj$seurat_clusters)]))
  temp_msi_score <- filter(cluster_info, clusters == temp_cluster)
   s_obj$sensor2[i] <- temp_msi_score$msi
  
}

return(s_obj)
  
  
}


add_pro <- function(s_obj,sample_name){
  s_obj[['pro']] <- NA
cluster_info <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'))
cluster_info <- data.frame(clusters = sort(unique(cluster_info$V2)))
        for(i in 1:nrow(cluster_info)){
   num <- cluster_info$clusters[i]
      temp15 <- fread(paste0('../pro_results/',sample_name,'_cluster_',num))
      cluster_info$msi[i] <- temp15$`%`
      cluster_info$tot_sites[i] <- temp15$Total_Number_of_Sites
        if(cluster_info$tot_sites[i] < 10 ){cluster_info$msi[i] <- NA}
}
for(i in 1:nrow(s_obj[['pro']])){
  temp_name <- names(s_obj$pro[i])
  temp_cluster <- as.numeric(as.character(s_obj$seurat_clusters[temp_name == names(s_obj$seurat_clusters)]))
  temp_msi_score <- filter(cluster_info, clusters == temp_cluster)
   s_obj$pro[i] <- temp_msi_score$msi

}

return(s_obj)


}

add_msings <- function(s_obj,sample_name){
s_obj[['msings']] <- NA
cluster_info <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'))
cluster_info <- data.frame(clusters = sort(unique(cluster_info$V2)))
for(i in 1:nrow(cluster_info)){
num <- cluster_info$clusters[i]
temp <- fread(paste0('../msings_results/',sample_name,'_cluster_',num,'.MSI_Analysis.txt'))
cluster_info$msi[i] <- temp[3,2]
cluster_info$tot_sites[i] <- temp[2,2]
if(cluster_info$tot_sites[i] < 10 ){cluster_info$msi[i] <- NA}
}
for(i in 1:nrow(s_obj[['msings']])){
  temp_name <- names(s_obj$msings[i])
  temp_cluster <- as.numeric(as.character(s_obj$seurat_clusters[temp_name == names(s_obj$seurat_clusters)]))
  temp_msi_score <- filter(cluster_info, clusters == temp_cluster)
   s_obj$msings[i] <- as.numeric(temp_msi_score$msi)*100

}

return(s_obj)


}



add_numbat <- function(s_obj, sample_name){
nb = Numbat$new(out_dir = paste0('../numbat/',sample_name))
clone_assignments <- data.frame(cell=nb$clone_post$cell, type=nb$clone_post$compartment_opt)

s_obj[['numbat']] <- NA

for(i in 1:nrow(s_obj[['numbat']])){
temp_name <- names(s_obj$numbat[i])
temp_clone <- filter(clone_assignments, cell == temp_name)

if(nrow(temp_clone)==0){s_obj$numbat[i] <- NA} else{s_obj$numbat[i] <- temp_clone$type}
}
return(s_obj)
}


add_cc <- function(s_obj,sample_name){
copykat <- fread(paste0('../copy_kat_results/',sample_name,'_copykat_prediction.txt'))
s_obj[['copykat']] <- NA 
for (i in 1:length(names(s_obj$copykat))){
temp2 <- filter(copykat, cell.names == names(s_obj$copykat[i]))
s_obj$copykat[i] <- temp2$copykat.pred
}
return(s_obj)
}

add_muts <- function(s_obj,sample_name){
cluster_info <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'))
cluster_info <- data.frame(clusters = sort(unique(cluster_info$V2)))
for(i in 1:nrow(cluster_info)){
num <- cluster_info$clusters[i]
temp <- readRDS(paste0('../vartrix_results/',sample_name,'.',num,'.rds'))
s_obj <- AddMetaData(object=s_obj,metadata=t(temp))
}
return(s_obj)

}

add_muts2 <- function(s_obj,sample_name){

cluster_info <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'))
colnames(cluster_info) <- c('barcode','cluster')

maf_files = list.files(path =paste0("../maf/",sample_name,"/"), pattern = ".maf", full.names = TRUE)
merged_maf <-maftools::merge_mafs(mafs = maf_files, verbose = TRUE) 
maf_mat <- matrix(data=NA,ncol=length(merged_maf@data$Hugo_Symbol),nrow=(length(colnames(s_obj))))
colnames(maf_mat) <- merged_maf@data$Hugo_Symbol
rownames(maf_mat) <- colnames(s_obj)

n_cores=detectCores()-1

cellular <- function(cells){


cell <- rownames(maf_mat)[cells]
temp_mat <- t(as.matrix(maf_mat[cell,],colnames=TRUE))
rownames(temp_mat) <- cell
temp_info <- filter(cluster_info, barcode == cell)
num <- temp_info$cluster
temp_maf <- read.maf(paste0('../maf/',sample_name,'/',sample_name,'.',num,'.maf'))


for(genes in 1:ncol(maf_mat)){
gene <- colnames(maf_mat)[genes]
tiny <- filter(temp_maf@data, Hugo_Symbol == gene)

if(nrow(tiny) > 1){

temp_mat[,gene] <- max(tiny$t_depth)

}
}

if(cells %% 5 == 0){
message(paste0('finished cell ', cell, '; ', cells/nrow(maf_mat)*100,'% complete'))
}


return(temp_mat)
}

matrix_list <- mclapply(1:nrow(maf_mat), cellular,mc.cores=n_cores)

maf_mat2 <- do.call(rbind,matrix_list[1:length(matrix_list)])
#filter out genes that only have NA values for each cell
maf_mat_trim <- maf_mat2[, colMaxs(maf_mat2,na.rm=TRUE) > 15]
#maf_mat_trim <- maf_mat2[,colSums(is.na(maf_mat))<nrow(maf_mat2)]

s_obj <- AddMetaData(object=s_obj,metadata=maf_mat_trim)
return(s_obj)
}



sample <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))


new_sample <- add_sensor2(sample,sample_name)
new_sample2 <- add_numbat(new_sample,sample_name)
new_sample3 <- add_cc(new_sample2,sample_name)
new_sample4 <- add_pro(new_sample3,sample_name)
new_sample5 <- add_msings(new_sample4,sample_name)
new_sample6 <- add_muts2(new_sample5,sample_name)


saveRDS(new_sample6, paste0('../annotated_h5/',sample_name,'.rds'))

