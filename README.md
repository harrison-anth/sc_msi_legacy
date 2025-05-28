# Intratumoral heterogeneity in microsatellite instability status at single cell resolution
## Article DOI: Under review

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie

### License type
MIT; see LICENSE file for more information.

### Repository information

This repository contains all the code needed to reproduce the results of our recent manuscript and will be archived after the manuscript is published.
Also included in this repository are the raw results (MSI scores, summary statistics, etc.).
To access the distributable version of the SC-MSI pipeline visit https://github.com/harrison-anth/sc_msi

Note: Each directory in this repository might have a deprecated code or files folder. These, while useful to archive, were not used in the final results for the manuscript.

## Information on directories and files in this repository

### conda_envs
Contains all conda environments and Snakemake profiles

atomic.yml -- used for any R code that requires scATOMIC to be run

seurat.yml -- used for any R code that involves Seurat objects

slurm_executor_profile -- contains settings for using the SLURM executor plugin with Snakemake

### images
Contains preliminary plots for all the samples prior to integration 

### manifests
This directory contains all the manifest files used to run the Snakemake pipeline.

all_artificial_samples_{1..10}.tsv -- the mix ID's for different mixing experiment runs.

all_gsm_samples.txt -- sample names for data from GSE205506

all_mixing_tables_manifest.tsv -- master manifest file to run the Snakemake pipeline for mixed samples

all_mix_patients.txt -- formatted mix ID's for Snakemake pipeline

all_samples.tsv -- All samples for data from EGAD00001008555, EGAD00001008584, EGAD00001008585, and PRJNA932556

combined_key.tsv -- A key of all samples and patient IDs for use with the Snakemake pipeline

excluded_sample_names.txt -- list of excluded sample names removed from the master key due to no identified cancer cells

final_key.tsv -- key that links sample names and patient id's

glob_patients.txt -- list of all patients (for use with MSIsensor-RNA glob script)

gsm_patients.txt -- list of patients from GSE205506

mixing_table_manifest_{1..100).tsv -- information on the proportion of cells sampled for each mixing run

mix_key{1..10}.tsv -- keys used to link file name and mix ID for Snakemake pipeline

mix_patients{1..10}.txt -- list of mix ID's for sanekamke pipeline

patient_ids.txt -- list of patients from EGAD00001008555, EGAD00001008584, EGAD00001008585, and PRJNA932556

### markdown_files
parallelized_reporter.rmd -- Patient report generator

### mix_summary_stats
The summary stats files for all mixes in the study

### msings_results
mSINGS scores for each individual in the study with a BAM file that was run at the individual and pseudobulk level (left out of final results due to lack of agreement with scATOMIC)

### pro_results
MSIsensor-pro scores for each individual in the study with a BAM file that was run at the individual and pseudobulk level (left out of final results due to lack of agreement with scATOMIC)

### pseudobulk_barcodes
The barcodes for each individual/mix in the manuscript used to generate pseudobulk data

### scripts
annotate_bamog.R -- used to annotate seurat objects with MSI score, number of subclones, etc. from BAM files

annotate_gsmog.R -- used to annotate seurat objects with MSI score, number of subclones, etc. from MTX files





Please feel free to reach out with any questions if this README has not answered your questions. 

