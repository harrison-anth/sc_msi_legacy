#first iteration based on function found here: https://bioinformatics-core-shared-training.github.io/SingleCell_RNASeq_May23/UnivCambridge_ScRnaSeqIntro_Base/Markdowns/101-seurat_part2.html
#also see the merge function from the seurat v4.3 archive https://satijalab.org/seurat/archive/v4.3/merge

#set seed (remove this for random sampling)
set.seed(9)

#load libs
library(Seurat)
library(subSeq)
library(data.table)
library(tidyverse)
library(R.utils)

#add command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
manifest <- as.character(argus[1])

#load data 
key <- fread(paste0('../manifests',manifest))




#set seed
set.seed(seed = 152727)





#define funcs

#s_obj should be a seurat object
#prop is the proportion of cells to sample from this object
cell_sampler <- function(s_obj,prop){

num_cells <- ncol(s_obj)
cells <- num_cells * prop



subbed <- s_obj[, sample(colnames(s_obj), size = cells, replace=F)]

return(subbed)
}


merged <- merge(MSIH_samp, y = MSS_samp, add.cell.ids = c("MSIH", "MSS"), project = "genesis")

