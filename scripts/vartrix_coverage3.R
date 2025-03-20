library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)
library(gridExtra)

#good reference https://rnabio.org/module-08-scrna/0008/06/01/Cancer_cell_identification/

#define print_mut func
print_mut <- function(s_obj, genos_matrix, alt_matrix, ref_matrix, frac_matrix, barcodes){

for(i in 1:ncol(genos_matrix)){
column <- colnames(genos_matrix)[i]
column_df <- data.frame(genos_matrix[,i])
row.names(column_df) <- barcodes$V1
colnames(column_df) <- column

gene_name <- colnames(add_anno(column_df,annotations))

#make genotypes readable
column_df[,column] <- str_replace(as.character(column_df[,column]), "0", "No Call")
column_df[,column] <- str_replace(as.character(column_df[,column]), "1", "ref/ref")
column_df[,column] <- str_replace(as.character(column_df[,column]), "2", "alt/alt")
column_df[,column] <- str_replace(as.character(column_df[,column]), "3", "alt/ref")

s_obj <- AddMetaData(object=s_obj,metadata=column_df)

s1 <- DimPlot(s_obj,group.by=c(paste0(column)))+ggtitle(paste(gene_name,'Genotypes'))

alt_temp <- data.frame(alt_matrix[,column])
row.names(alt_temp) <- barcodes$V1
colnames(alt_temp) <- column

s_obj <- AddMetaData(object=s_obj,metadata=alt_temp)
s2 <- FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'alt counts'))

ref_temp <- data.frame(ref_matrix[,column])
row.names(ref_temp) <- barcodes$V1
colnames(ref_temp) <- column
s_obj <- AddMetaData(object=s_obj,metadata=ref_temp)
s3 <-  FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'ref counts'))

frac_temp <- data.frame(frac_matrix[,column])
row.names(frac_temp) <- barcodes$V1
colnames(frac_temp) <- column
s_obj <- AddMetaData(object=s_obj,metadata=frac_temp)
s4 <-FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'AF'))


print(marrangeGrob(list(s1,s2,s3,s3),nrow=2,ncol=2))

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

#add CCF function to estimate somatic mutation in cancer subpopulation

add_ccf <- function(){
set.seed(9)
n <- 1000
X <- runif(n, -2, 2)  # Predictor variable
beta_0 <- 0.5
beta_1 <- 2
c_val <- 0.8  # Upper asymptote (less than 1)
# Generate probabilities using scaled logistic function
P <- c_val / (1 + exp(-(beta_0 + beta_1 * X)))
# Generate binary outcomes
Y <- rbinom(n, 1, P)
# Fit the scaled logistic model using nls()
scaled_logistic_model <- nls(Y ~ c_param / (1 + exp(-(b0 + b1 * X))),
                             start = list(c_param = 0.8, b0 = 0, b1 = 1),
                             algorithm = "port",
                             lower = c(0, -Inf, -Inf),  # Ensure c is positive
                             upper = c(1, Inf, Inf))    # Ensure c is at most 1
# View results
summary(scaled_logistic_model)
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
#snv_matrix_alt <- snv_matrix_alt[colSums(snv_matrix_alt)>=nrow(snv_matrix_alt)]


snv_matrix_frac <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_frac_matrix.mmx'))))
colnames(snv_matrix_frac) <- barcodes$V1
row.names(snv_matrix_frac) <- snps$V1
snv_matrix_frac <- as.data.frame(t(snv_matrix_frac))

#while redundant to filter, it saves time.
snv_matrix_frac <- snv_matrix_frac[colSums(snv_matrix_frac,na.rm=TRUE)>=5]

snv_matrix_genos <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_genos_matrix.mmx'))))
colnames(snv_matrix_genos) <- barcodes$V1
row.names(snv_matrix_genos) <- snps$V1
snv_matrix_genos <- as.data.frame(t(snv_matrix_genos))
snv_matrix_genos <- snv_matrix_genos[colSums(snv_matrix_genos)>=nrow(snv_matrix_genos)]






####try rerunning vartrix with allele-fraction mode


#read in annotated seurat object
s_obj <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

#read in annovar annotations
annotations <- fread(paste0('../annovar_results/',sample_name,'.avinput.variant_function'))
annotations$key <- paste0(annotations$V3,':',annotations$V4)





#print it out all -- consider filtering further if too many snps or uninformative images are retained --
pdf(paste0('../images/',sample_name,'_muts.pdf'),onefile=TRUE,height=12,width=26)
print_mut(s_obj=s_obj,genos_matrix=snv_matrix_genos,alt_matrix=snv_matrix_alt,ref_matrix=snv_matrix_ref,frac_matrix=snv_matrix_frac,barcodes=barcodes)
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


