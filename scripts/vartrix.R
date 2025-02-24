library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

snps <- read.table(paste0("../vartrix_results/",sample_name,".snv_loci.txt"), header = F)
barcodes <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'.tsv'),header=FALSE)

snv_matrix <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'.mmx'))))
colnames(snv_matrix) <- barcodes$V1
row.names(snv_matrix) <- snps$V1



saveRDS(mat2, paste0('../vartrix_results/',sample_name,'.',num,'.rds'))


