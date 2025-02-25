library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

snps <- read.table(paste0("../vartrix_results/",sample_name,".snv_loci.txt"), header = F)
barcodes <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_just_barcodes.tsv'),header=FALSE)

snv_matrix <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'.mmx'))))
colnames(snv_matrix) <- barcodes$V1
row.names(snv_matrix) <- snps$V1


annotations <- fread(paste0('../annovar_results/',sample_name,'.',num,'.avinput.variant_function'))
annotations$key <- paste0(annotations$V3,':',annotations$V4)

for(i in 1:nrow(snv_matrix)){
temp2 <- filter(annotations, annotations$key == rownames(snv_matrix[i,]))

if(rownames(snv_matrix)[i] %in% annotations$key){

#prevent duplicate rownames
if(temp2$V2[1] %in% rownames(snv_matrix)){
rownames(snv_matrix)[i] <- paste0(temp2$V2[1],'_',i)
} else{

rownames(snv_matrix)[i] <- temp2$V2[1]}}

}

#list of cool mutations

muts <- c('APC','TP53','KRAS','BRAF','EGFR','MLH','MSH','PSM')
mat2 <- rbind(snv_matrix[grep(pattern='APC',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='TP53',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='KRAS',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='BRAF',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='EGFR',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='MLH',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='MSH',x=rownames(snv_matrix),ignore.case=TRUE),],
	snv_matrix[grep(pattern='PSM',x=rownames(snv_matrix),ignore.case=TRUE),])

saveRDS(mat2, paste0('../vartrix_results/',sample_name,'.',num,'.rds'))


