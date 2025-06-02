#annotate filtered h5 files
#functions and libs

library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)
library(data.table)
library(numbat)


argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

set.seed(seed = 152727)

sample <- Read10X_h5(paste0('../cell_ranger_output/',sample_name,'/outs/filtered_feature_bc_matrix.h5'))
allele_counts <- fread(paste0('gunzip -cq ../numbat/',sample_name,'/',sample_name,'_allele_counts.tsv.gz'))
out = run_numbat(sample,ref_hca,allele_counts,genome="hg38",t=1e-5,ncores=1,plot=TRUE,out_dir=paste0('../numbat/',sample_name,'/'))

