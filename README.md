# Intratumoral heterogeneity in microsatellite instability status at single cell resolution
## Article DOI: Under review

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie

### License type
MIT with attribution and commercial use restriction; see LICENSE file for more information.

### Repository information

This repository contains all the code needed to reproduce the results of our recent manuscript and will be archived after the manuscript is published.
Also included in this repository are the raw results (MSI scores, summary statistics, etc.).
To access the distributable version of the SC-MSI pipeline visit https://github.com/harrison-anth/sc_msi

## Information on directories and files in this repository
### Conda environments and Snakemake profiles
The conda_envs folder contains all conda environments and Snakemake profiles

atomic.yml -- used for any R code that requires scATOMIC to be run
seurat.yml -- used for any R code that involves Seurat objects
slurm_executor_profile -- contains settings for using the SLURM executor plugin with Snakemake

### images
Contains preliminary plots for all the samples prior to integration 

### manifests
This directory contains all the manifest files used to run the Snakemake pipeline
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





Please feel free to reach out with any questions if this README has not answered your questions. 

