library(data.table)
library(Matrix)
library(Seurat)
library(infercnv)
library(tidyverse)
library(R.utils)
library(stringr)

set.seed(seed = 152727)

argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
outdir <- as.character(argus[2])

s_obj <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))
#create annotations file
anno <- data.frame(cell_name=colnames(s_obj),cancer=s_obj$pan_cancer_cluster,
type=s_obj$scATOMIC_pred,msi=s_obj$sensor_rna_status)

#change MSI NA's to MSS (on the basis that MSS unless otherwise proven high)
anno$msi[is.na(anno$msi)] <- "MSS"

anno$cancer2 <- paste0(anno$cancer,'_',anno$msi)
final_anno <- data.frame(cell_name=anno$cell_name,type=ifelse(anno$cancer =="Normal",yes=anno$cancer,no=anno$cancer2))
fwrite(final_anno,file=paste0('../temp/',sample_name,'_anno.tsv'),sep='\t',col.names=FALSE)

#create infercnv object
cnv_obj <- CreateInfercnvObject(raw_counts_matrix=s_obj@assays$RNA$counts,
annotations_file=paste0('../temp/',sample_name,'_anno.tsv'),
delim='\t',gene_order_file='../temp/hg38_gencode_v27.txt',
ref_group_names=c("Normal"))

#fix scipen for infercnv run
options(scipen = 100)

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(cnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(outdir),
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
                             )
