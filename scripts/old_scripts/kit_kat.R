library(Seurat) 
library(patchwork) 
library(scAnnotatR) 
library(PreMSIm)
library(copykat)
library(R.utils)
library(data.table)
library(tidyverse)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
core_num <- as.numeric(argus[2])
gsm <- as.character(argus[3])
gsm_key <- fread('../manifests/final_gsm_key.tsv')


set.seed(seed = 152727)

if(gsm == 'Y'){

indv_key <- filter(gsm_key, sample_id == sample_name)

temp <- ReadMtx(
  mtx=paste0('../gsm_samps/',indv_key$file_prefix,'_matrix.mtx.gz'),
  cells=paste0('../gsm_samps/',indv_key$file_prefix,'_barcodes.tsv.gz'),
  features=paste0('../gsm_samps/',indv_key$file_prefix,'_features.tsv.gz')
)} else{
temp <- Read10X_h5(paste0('../cell_ranger_output/',sample_name,'/outs/filtered_feature_bc_matrix.h5'))
}

temp <-  CreateSeuratObject(counts=temp,
                          project="MSI",
                          min.cells=0,
                          min.features=0)
raw_data<- GetAssayData(temp,assay='RNA',layer='counts')

copykat.test <- copykat(rawmat=raw_data, id.type="S",
ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=sample_name, 
distance="euclidean", norm.cell.names="",output.seg="FLASE", 
plot.genes="FALSE", 
genome="hg20",
n.cores=core_num)
