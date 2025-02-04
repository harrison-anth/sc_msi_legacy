#!/bin/bash

#mix proportions of cells from different single-cell RNA sequencing samples together in R

#declare variables

#key with sample names and PCR/IHC ground truth. 


manifest=../manifests/final_key.tsv


module load Anaconda3
conda activate seurat

Rscript cellMixe.R $manifest



