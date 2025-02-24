#!/bin/bash

#mix proportions of cells from different single-cell RNA sequencing samples together in R

#declare variables

#key with sample names and PCR/IHC ground truth. 


manifest=../manifests/final_key.tsv


module load Anaconda3
conda activate seurat

#msih homogenous
sample1_name='GSM6213995'
#mss homogenousish
sample2_name='XHC118-SI-GA-F1'

Rscript cellMixe.R $manifest $sample1_name $sample2_name



