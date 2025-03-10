library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)

#good reference https://rnabio.org/module-08-scrna/0008/06/01/Cancer_cell_identification/

#define print_mut func
print_mut <- function(s_obj, snv_matrix,barcodes){
for(i in 1:ncol(snv_matrix)){
column <- colnames(snv_matrix)[i]
column_df <- data.frame(snv_matrix[,i])
row.names(column_df) <- barcodes$V1
colnames(column_df) <- column
s_obj <- AddMetaData(object=s_obj,metadata=column_df)
 print(DimPlot(s_obj,group.by=c(paste0(column))))
}
}




argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

snps <- read.table(paste0("../vartrix_results/",sample_name,".snv_loci.txt"), header = F)
barcodes <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/','just_barcodes.tsv'),header=FALSE)

snv_matrix <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_alt_matrix.mmx'))))
colnames(snv_matrix) <- barcodes$V1
row.names(snv_matrix) <- snps$V1

snv_matrix <- as.data.frame(t(snv_matrix))

s_obj <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))

#filter out sites with too much missing data. Simple filter. Because 0 is a no-call
#Good to have at least 1 call on average per cell, but the 2 and 3 genos will allow for some.
snv_filt <- snv_matrix[colSums(snv_matrix)>=nrow(snv_matrix)]

#print it out all -- consider filtering further if too many snps or uninformative images are retained --
pdf(paste0('../images/',sample_name,'_muts.pdf'))
print_mut(s_obj=s_obj,snv_matrix=snv_filt,barcodes=barcodes)
dev.off()



#optional annotation code; not yet implemented

#annotations <- fread(paste0('../annovar_results/',sample_name,'.',num,'.avinput.variant_function'))
#annotations$key <- paste0(annotations$V3,':',annotations$V4)

#for(i in 1:nrow(snv_matrix)){
#temp2 <- filter(annotations, annotations$key == rownames(snv_matrix[i,]))

#if(rownames(snv_matrix)[i] %in% annotations$key){

#prevent duplicate rownames
#if(temp2$V2[1] %in% rownames(snv_matrix)){
#rownames(snv_matrix)[i] <- paste0(temp2$V2[1],'_',i)
#} else{
#
#rownames(snv_matrix)[i] <- temp2$V2[1]}}
#
#}
#
#list of cool mutations

#muts <- c('APC','TP53','KRAS','BRAF','EGFR','MLH','MSH','PSM')
#mat2 <- rbind(snv_matrix[grep(pattern='APC',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='TP53',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='KRAS',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='BRAF',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='EGFR',x=rownames(snv_matrix),ignore.case=TRUE),],
##	snv_matrix[grep(pattern='MLH',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='MSH',x=rownames(snv_matrix),ignore.case=TRUE),],
#	snv_matrix[grep(pattern='PSM',x=rownames(snv_matrix),ignore.case=TRUE),])
#
#saveRDS(mat2, paste0('../vartrix_results/',sample_name,'.',num,'.rds'))


