#setwd
#lobal libs
library(Seurat)
library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(R.utils)
library(gridExtra)

#create not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

#get sample name
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])




#define funk
cnv_scannerprint_mut <- function(s_obj){

    #add ccf information for each cluster
    for(L in levels(s_obj$seurat_clusters)){
    mini_s_obj <- subset(s_obj, seurat_clusters == L)

    mini_binary <- as.integer(as.logical(s_obj$has_loss_chr1))
    mini_depth <- rowSums(t(mini@assays$RNA$counts)))

     ccf_model_summary <- summary(ccf_model)
      column_df[rownames(column_df)%in% colnames(mini_s_obj),"ccf"] <- ccf_model_summary$coefficients[1]

      column_df[rownames(column_df)%in% colnames(mini_s_obj),"p_value"] <- ccf_model_summary$coefficients[10]
      
    }}

#add CCF function to estimate somatic mutation in cancer subpopulation

add_ccf <- function(depth_df,binary_df){
  set.seed(9)
  X <- depth_df
  beta_0 <- 0.5
  beta_1 <- 2
  c_val <- 0.8  # Upper asymptote (less than 1)
  # Generate probabilities using scaled logistic function
  P <- c_val / (1 + exp(-(beta_0 + beta_1 * X)))
  # Generate binary outcomes
  Y <- as.numeric(binary_df) #whether snp was detected; bin
  # Fit the scaled logistic model using nls()
  scaled_logistic_model <-  nls(Y ~ c_param / (1 + exp(-(b0 + b1 * X))),
                                              start = list(c_param = 0.8, b0 = 0, b1 = 1),
                                              algorithm = "port",
                                              lower = c(0, -Inf, -Inf),  # Ensure c is positive
                                              upper = c(1, Inf, Inf))
                                     
  
  return((scaled_logistic_model))}










snps <- read.table(paste0('../vartrix_results/',sample_name,".snv_loci.txt"), header = F)
barcodes <- fread(paste0('../pseudobulk_barcodes/',sample_name,'/just_barcodes.tsv'),header=FALSE)

snv_matrix_genos <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_genos_matrix.mmx'))))
colnames(snv_matrix_genos) <- barcodes$V1
row.names(snv_matrix_genos) <- snps$V1
snv_matrix_genos <- snv_matrix_genos[rowSums(snv_matrix_genos)>=50,]
#gc()
snv_matrix_genos <- as.data.frame(t(snv_matrix_genos))


snv_matrix_alt <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_alt_matrix.mmx'))))
colnames(snv_matrix_alt) <- barcodes$V1
row.names(snv_matrix_alt) <- snps$V1
snv_matrix_alt <- snv_matrix_alt[rowSums(snv_matrix_alt)>=100,]
#gc()
snv_matrix_alt <- as.data.frame(t(snv_matrix_alt))


snv_matrix_ref <- as.data.frame(as.matrix(readMM(paste0('../vartrix_results/',sample_name,'_ref_matrix.mmx'))))
colnames(snv_matrix_ref) <- barcodes$V1
row.names(snv_matrix_ref) <- snps$V1
snv_matrix_ref <- snv_matrix_ref[rowSums(snv_matrix_ref)>=1,]
snv_matrix_ref <- as.data.frame(t(snv_matrix_ref))
#gc()

#snv_matrix_frac <- as.data.frame(as.matrix(readMM(paste0(sample_name,'_frac_matrix.mmx'))))
#colnames(snv_matrix_frac) <- barcodes$V1
#row.names(snv_matrix_frac) <- snps$V1
#snv_matrix_frac <- snv_matrix_frac[rowSums(snv_matrix_frac)>=10,]
##snv_matrix_genos <- as.data.frame(t(snv_matrix_genos))
#gc()

#ensure commonalities

alt_matrix <- snv_matrix_alt
#[colnames(snv_matrix_alt)%in% colnames(snv_matrix_ref)]
ref_matrix <- snv_matrix_ref
#[colnames(snv_matrix_ref)%in% colnames(snv_matrix_alt)]
genos_matrix <- snv_matrix_genos[colnames(snv_matrix_genos)%in% colnames(snv_matrix_alt)]








s_obj <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

#read in annovar annotations
annotations <- fread(paste0('../annovar_results/' ,sample_name,'.avinput.variant_function'))
annotations$key <- paste0(annotations$V3,':',annotations$V4)

pdf(paste0('../images/',sample_name,'_ccf.pdf'),onefile=TRUE,height=12,width=26)
print_mut(s_obj=s_obj,genos_matrix=genos_matrix,alt_matrix=alt_matrix,ref_matrix=ref_matrix,barcodes=barcodes)
dev.off()






