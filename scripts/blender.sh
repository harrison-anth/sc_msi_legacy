#!/bin/bash

#mix proportions of cells from different single-cell RNA sequencing samples together in R

#declare variables

#key with sample names and PCR/IHC ground truth. 


manifest=../manifests/final_key.tsv


module load Anaconda3
conda activate seurat

sample1_name='GSM6213971'
sample2_name='GSM6213963'

Rscript cellMixe.R $manifest $sample1_name $sample2_name



