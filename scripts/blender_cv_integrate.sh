#!/bin/bash

#mix proportions of cells from different single-cell RNA sequencing samples together in R

#number of cross validations (runs)
cv=10

#key with sample names and PCR/IHC ground truth. 
manifest=../manifests/final_key3.tsv


module load Anaconda3
conda activate seurat

#msih homogenous
sample1_name='GSM6213995'
#mss homogenousish
sample2_name='XHC118-SI-GA-F1'

Rscript cellMixeR_cv_integrate.R $manifest $sample1_name $sample2_name $cv



