library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)

#good reference https://rnabio.org/module-08-scrna/0008/06/01/Cancer_cell_identification/

#define print_mut func
print_mut <- function(s_obj, alt_matrix, ref_matrix, frac_matrix, barcodes){

for(i in 1:ncol(alt_matrix)){
column <- colnames(alt_matrix)[i]
column_df <- data.frame(alt_matrix[,i])
row.names(column_df) <- barcodes$V1
colnames(column_df) <- column

gene_name <- colnames(add_anno(column_df,annotations))

s_obj <- AddMetaData(object=s_obj,metadata=column_df)
 print(FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'alt counts')))

ref_temp <- data.frame(ref_matrix[,column])
row.names(ref_temp) <- barcodes$V1
colnames(ref_temp) <- column
s_obj <- AddMetaData(object=s_obj,metadata=ref_temp)
 print(FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'ref counts')))

frac_temp <- data.frame(frac_matrix[,column])
row.names(frac_temp) <- barcodes$V1
colnames(frac_temp) <- column
s_obj <- AddMetaData(object=s_obj,metadata=frac_temp)
 print(FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'AF')))

}
}

#add function to add annotations to snv_matrix
add_anno <- function(snv_matrix,annotations){

for(i in 1:ncol(snv_matrix)){
temp2 <- filter(annotations, annotations$key == colnames(snv_matrix)[i])

if(colnames(snv_matrix)[i] %in% annotations$key){

#prevent duplicate rownames
if(temp2$V2[1] %in% colnames(snv_matrix)){
colnames(snv_matrix)[i] <- paste0(temp2$V2[1],'_',i)
} else{

colnames(snv_matrix)[i] <- temp2$V2[1]}}

}
return(snv_matrix)
}


argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

snps <- read.table(paste0("../vartrix_results/",sample_name,".snv_loci.txt"), header = F)
barcodes <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/','just_barcodes.tsv'),header=FALSE)

snv_matrix_alt <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_alt_matrix.mmx'))))

colnames(snv_matrix_alt) <- barcodes$V1
row.names(snv_matrix_alt) <- snps$V1
snv_matrix_alt <- as.data.frame(t(snv_matrix_alt))

snv_matrix_ref <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_ref_matrix.mmx'))))
colnames(snv_matrix_ref) <- barcodes$V1
row.names(snv_matrix_ref) <- snps$V1
snv_matrix_ref <- as.data.frame(t(snv_matrix_ref))


#hold off on filtering ref. Might have homozygous for alt
#snv_matrix_ref <- snv_matrix_ref[colSums(snv_matrix_ref)>=nrow(snv_matrix_ref)]

#filter based on alt depth
snv_matrix_alt <- snv_matrix_alt[colSums(snv_matrix_alt)>=nrow(snv_matrix_alt)]


snv_matrix_frac <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_frac_matrix.mmx'))))
colnames(snv_matrix_frac) <- barcodes$V1
row.names(snv_matrix_frac) <- snps$V1
snv_matrix_frac <- as.data.frame(t(snv_matrix_frac))

#while redundant to filter, it saves time.
snv_matrix_frac <- snv_matrix_frac[colSums(snv_matrix_frac,na.rm=TRUE)>=5]



####try rerunning vartrix with allele-fraction mode


#read in annotated seurat object
s_obj <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

#read in annovar annotations
annotations <- fread(paste0('../annovar_results/',sample_name,'.avinput.variant_function'))
annotations$key <- paste0(annotations$V3,':',annotations$V4)






#print it out all -- consider filtering further if too many snps or uninformative images are retained --
pdf(paste0('../images/',sample_name,'_muts.pdf'))
print_mut(s_obj=s_obj,alt_matrix=snv_matrix_alt,ref_matrix=snv_matrix_ref,frac_matrix=snv_matrix_frac,barcodes=barcodes)
dev.off()








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


