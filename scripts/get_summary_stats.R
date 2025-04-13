#Integrate sample from common patient ID's
#Much of the integration code from user alg697's comment on seurat github issue
#https://github.com/satijalab/seurat/issues/8938

#functions and libs
library(data.table)
library(Matrix)
library(Seurat)
library(R.utils)
library(tidyverse)
library(infercnv)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
gsm <- as.character(argus[2])

set.seed(seed = 152727)

if(gsm == 'Y'){key <- fread('../manifests/final_gsm_key.tsv',header=TRUE)
} else if(gsm == 'N'){key <- fread('../manifests/final_key3.tsv',header=TRUE)}

get_summary_stats <- function(sample_name){

int_s_obj <-readRDS(paste0('../integrated_samples/',sample_name,'.rds'))
int_s_obj <- add_to_seurat(seurat_obj=int_s_obj, 
infercnv_output_path=paste0('/data3/hanthony/infer_cnv_temp/patient_',sample_name,'/'),top_n=10)

int_s_obj <- subset(int_s_obj, subset=tissue != "Normal" & tissue != "normal")

int_cancer_obj <- readRDS(paste0('../integrated_samples/',sample_name,'_cancer.rds'))

cancer_fin <- AddMetaData(int_cancer_obj,metadata=int_s_obj$infercnv_subcluster,col.name='infercnv_subcluster')

indv_key <- filter(key, patient_id == sample_name)
num_tumor_samps <- nrow(filter(indv_key, site != "normal" & site != "Normal"))
num_tot_samps <- nrow(indv_key)


if(length(levels(droplevels(cancer_fin$seurat_clusters)))<2){
print('Only one cluster of cancer cells; cannot run anova')
clustys <- cancer_fin
anova_df <- data.frame(patient_id=sample_name,
			DF=0,
			Ssq="NA",
			Msq="NA",
			F="NA",
			P="NA")
subclone_df <- data.frame(Patient=sample_name,
			Clusters=1,
			F=NA,
			NumCc = ncol(cancer_fin),
			Num_Samp=num_tumor_samps,
                        Num_MSS_subclones=length(na.omit(str_extract(string=unique(cancer_fin$infercnv_subcluster),pattern='MSS'))),
                        Num_MSIH_subclones=length(na.omit(str_extract(string=unique(cancer_fin$infercnv_subcluster),pattern='MSI-H'))),
                        Tot_subclones=length(unique(cancer_fin$infercnv_subcluster)))

fwrite(x=anova_df,file=paste0('../summary_stats/',sample_name,'_anova_results.tsv'),sep='\t')
fwrite(x=subclone_df,file=paste0('../summary_stats/',sample_name,'_cluster_stats.tsv'),sep='\t')





} else{

#Write out cluster f and summary stats
anova_df <- data.frame(cluster=as.character(cancer_fin$seurat_clusters),
                        msi_score=as.numeric(cancer_fin$sensor_rna_prob))
anova_result <- aov(msi_score ~ cluster,data=anova_df)

anova_df <- data.frame(patient_id=sample_name,
                        Df=summary(anova_result)[[1]][["Df"]][1],
                        Ssq=round(summary(anova_result)[[1]][["Sum Sq"]][1],4),
                        Msq=round(summary(anova_result)[[1]][["Mean Sq"]][1],4),
                        F = round(summary(anova_result)[[1]][["F value"]][1],4),
                        P = round(summary(anova_result)[[1]][["Pr(>F)"]][1],4))

subclone_df <- data.frame(Patient=sample_name,
                        Clusters=summary(anova_result)[[1]][["Df"]][1]+1,
                        F = round(summary(anova_result)[[1]][["F value"]][1],4),
                        numCc = ncol(cancer_fin),
                        Num_Samp=num_tumor_samps,
                        Num_MSS_subclones=length(na.omit(str_extract(string=unique(cancer_fin$infercnv_subcluster),pattern='MSS'))),
                        Num_MSIH_subclones=length(na.omit(str_extract(string=unique(cancer_fin$infercnv_subcluster),pattern='MSI-H'))),
                        Tot_subclones=length(unique(cancer_fin$infercnv_subcluster)))



fwrite(x=anova_df,file=paste0('../summary_stats/',sample_name,'_anova_results.tsv'),sep='\t')
fwrite(x=subclone_df,file=paste0('../summary_stats/',sample_name,'_cluster_stats.tsv'),sep='\t')




}
}





get_summary_stats(sample_name)

