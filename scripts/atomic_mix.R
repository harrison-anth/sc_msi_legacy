#annotate filtered h5 files
#load libraries
library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)
library(data.table)

#add command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

#set seed
set.seed(seed = 152727)

#Some issues with using ATOMIC on mix data (it doesn't seem able to label cells as cancerous if mixed between the two)
#using slightly dirty method, but it makes sense as the MSIsensor scores do not change from mix to mix. I think it's okay
#to keep the cell type annotations from the original samples as well. 


mix <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

atomic_cells2 <- readRDS(paste0('../atomic/','XHC118-SI-GA-F1','.rds'))
atomic_cells1 <- readRDS(paste0('../atomic/','GSM6213995','.rds'))
mergy <- merge(atomic_cells1,atomic_cells2)
mix$pan_cancer_cluster <- mergy$pan_cancer_cluster
mix$classification_confidence <- mergy$classification_confidence
mix$scATOMIC_pred <- mergy$scATOMIC_pred



saveRDS(mix, paste0('../atomic/',sample_name,'.rds'))

barcodes <- data.frame(barcode=names(mix$pan_cancer_cluster),call=mix$pan_cancer_cluster)

fwrite(barcodes,sep='\t',paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_cancer_barcodes.tsv'),
row.names=FALSE,col.names=FALSE)





