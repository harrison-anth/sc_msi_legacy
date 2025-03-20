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






#define print_mut func
print_mut <- function(s_obj, genos_matrix, alt_matrix, ref_matrix,barcodes){
  
  for(i in 1:ncol(genos_matrix)){
    column <- colnames(genos_matrix)[i]
    column_df <- data.frame(genos_matrix[,i])
    row.names(column_df) <- barcodes$V1
    colnames(column_df) <- column
    gene_name <- colnames(add_anno(column_df,annotations))
    
    alt_temp <- data.frame(alt_matrix[,column])
    row.names(alt_temp) <- barcodes$V1
    colnames(alt_temp) <- column
    
    s_obj <- AddMetaData(object=s_obj,metadata=alt_temp)
    s2 <- FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'alt counts'))
    
    
    if(column %!in% colnames(ref_matrix)){
      ref_temp <- alt_temp 
      ref_temp[,column] <- 0
      }else{    ref_temp <- data.frame(ref_matrix[,column])
    row.names(ref_temp) <- barcodes$V1
    colnames(ref_temp) <- column
      }
    s_obj <- AddMetaData(object=s_obj,metadata=ref_temp)
    s3 <-  FeaturePlot(s_obj,features=c(paste0(column)))+ggtitle(paste(gene_name,'ref counts'))
    
    column_df_binary <- data.frame(bin=ifelse(test = (column_df[,column]=="No Call" | column_df[,column]== "ref/ref"),yes="0",no="1"))    
    row.names(column_df_binary) <- rownames(ref_temp)

    column_df_depth <- alt_temp+ref_temp
    column_df$ccf <-NA
    column_df$p_value <- NA

#correct genotypes based off depth (i.e. can't genotype without appropriate reads)
for(C in 1:nrow(column_df)){
if(column_df_depth[C,column] <10){
column_df[C,column] <- 0}
}


#make genotypes readable
    column_df[,column] <- str_replace(as.character(column_df[,column]), "0", "No Call")
    column_df[,column] <- str_replace(as.character(column_df[,column]), "1", "ref/ref")
    column_df[,column] <- str_replace(as.character(column_df[,column]), "2", "alt/alt")
    column_df[,column] <- str_replace(as.character(column_df[,column]), "3", "alt/ref")

    s_obj <- AddMetaData(object=s_obj,metadata=column_df)

    s1 <- DimPlot(s_obj,group.by=c(paste0(column)))+ggtitle(paste(gene_name,'Genotypes'))








    #add ccf information for each cluster
    for(L in levels(s_obj$seurat_clusters)){
    mini_s_obj <- subset(s_obj, seurat_clusters == L)
    mini_binary <- as.numeric(column_df_binary[rownames(column_df_binary) %in% colnames(mini_s_obj),])
    mini_depth <- column_df_depth[rownames(column_df_depth) %in% colnames(mini_s_obj),]
    
    if(sum(mini_depth) < 10 | sum(mini_binary) <15){
      column_df[rownames(column_df)%in% colnames(mini_s_obj),"ccf"] <- 0
      column_df[rownames(column_df)%in% colnames(mini_s_obj),"p_value"] <- 0
      
      }else{ ccf_model <- add_ccf(depth_df =mini_depth, binary_df=mini_binary)
      if(sum(ccf_model[1] == "NA") == 10){
        column_df[rownames(column_df)%in% colnames(mini_s_obj),"ccf"] <- NA
        
        column_df[rownames(column_df)%in% colnames(mini_s_obj),"p_value"] <- NA
        
      }else if(ccf_model$message == "absolute function convergence (6)"){
        column_df[rownames(column_df)%in% colnames(mini_s_obj),"ccf"] <- NA
        
        column_df[rownames(column_df)%in% colnames(mini_s_obj),"p_value"] <- NA
      }
      else{
        ccf_model_summary <- summary(ccf_model)
      column_df[rownames(column_df)%in% colnames(mini_s_obj),"ccf"] <- ccf_model_summary$coefficients[1]

      column_df[rownames(column_df)%in% colnames(mini_s_obj),"p_value"] <- ccf_model_summary$coefficients[10]
      
    }}
    
    
    
    
    }
if(all(is.na(column_df$ccf))){
    print(marrangeGrob(list(s1,s2,s3,s3),nrow=2,ncol=2))
} else{

    s_obj <- AddMetaData(object=s_obj,metadata=column_df$ccf,col.name=paste0(column,'_ccf'))    
    s_obj <- AddMetaData(object=s_obj,metadata=column_df$p_value,col.name=paste0(column,'_pval'))    
    s4 <-  FeaturePlot(s_obj,features=c(paste0(column,'_ccf')))+ggtitle(paste(gene_name,' ccf'))
    s5 <- FeaturePlot(s_obj,features=c(paste0(column,'_pval')))+ggtitle(paste(gene_name,' pvals'))
    
    print(marrangeGrob(list(s1,s2,s3,s3,s4,s5),nrow=2,ncol=2))
}
    
  }
}





#add function to add annotations to snv_matrix
add_anno <- function(snv_matrix,annotations){
  
  for(A in 1:ncol(snv_matrix)){
    temp2 <- filter(annotations, annotations$key == colnames(snv_matrix)[A])
    
    if(colnames(snv_matrix)[A] %in% annotations$key){
      
      #prevent duplicate rownames
      if(temp2$V2[1] %in% colnames(snv_matrix)){
        colnames(snv_matrix)[A] <- paste0(temp2$V2[1],'_',A)
      } else{
        
        colnames(snv_matrix)[A] <- temp2$V2[1]}}
    
  }
  return(snv_matrix)
}




#add CCF function to estimate somatic mutation in cancer subpopulation

add_ccf <- function(depth_df,binary_df){
  set.seed(9)
  X <- depth_df
  #[,1]  # Predictor variable (depth)
  beta_0 <- 0.5
  beta_1 <- 2
  c_val <- 0.8  # Upper asymptote (less than 1)
  # Generate probabilities using scaled logistic function
  P <- c_val / (1 + exp(-(beta_0 + beta_1 * X)))
  # Generate binary outcomes
  Y <- as.numeric(binary_df) #whether snp was detected; bin
  # Fit the scaled logistic model using nls()
  scaled_logistic_model <-  tryCatch(expr=nls(Y ~ c_param / (1 + exp(-(b0 + b1 * X))),
                                              start = list(c_param = 0.8, b0 = 0, b1 = 1),
                                              algorithm = "port",
                                              lower = c(0, -Inf, -Inf),  # Ensure c is positive
                                              upper = c(1, Inf, Inf)),error=function(e)return(data.frame(coefficients=rep("NA",10))))
                                     
  
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






