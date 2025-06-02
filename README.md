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

slurm_profile -- contains the settings used specifically for the University of Galway HPC

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

artificial_config.yaml -- The config file used to run the Snakemake pipeline with artificial mix data

atomic_mix.R -- R script used to run scATOMIC on the artificial mix data

atomic.R -- R script used to run scATOMIC on BAM/MTX data

atomic.sh -- bash script used to call scATOMIC

bam_config.yaml -- configuration file for Snakemake used to handle BAM/FASTQ files

barcode_generator.R -- R script used to create and filter sample level seurat objects. Also exports cell barcode information to pseudobulk_barcodes/ directory

blender_cv.sh -- Bash script used to create artificial mix data

cellMix.R -- R script used to mix together different single-cell sequencing samples. (Currently hardcoded to integrate the expression data of each sample, but will make optional in the future).

get_summary_stats_mix.R -- R script to retrieve summary stats for mixed samples (number of subclone, MSI-H cells, etc.)

get_summary_stats.R -- R script to retrieve summary stats for a patient (number of subclones, MSI-H cells, etc.); mainly used in the final patient report.

glob_sensor_cancer.sh -- Quick bash script used to get the aggregated gene expression MSIsensor-RNA result for only the cancer cells in an individual

glob_sensor.sh -- Quick bash script used to get the aggregated gene expression MSIsensor-RNA result for all cells in an individual

glob_shaper_cancer.R -- Rscript used to arrange data into MSIsensor-RNA compatible structure (for use with an individual's aggregated cancer cells)

glob_shaper.R -- Rscript used to arrange data into MSIsensor-RNA compatible structure (for use with an individual's aggregated cells)

handle_bams.snake -- Snakefile used to run the SC-MSI pipeline on FASTQ/BAM files

handle_mix.snake -- Snakefile used to run the SC-MSI pipeline on mixture files 

infer_cnv.R -- R script used to run inferCNV on samples

infer_cnv_CRC2821.R -- R script used for CRC2821 inferCNV run (required a higher k_nn setting)

infer_patient_cnv.R -- R script used to run inferCNV on individuals

integrate_all.R -- R script that handles the integration step of the Seurat pipeline

integrate_mixes.R -- R script that handles the integration step of the Seurat pipeline for artificial mix samples

integrate.sh -- Bash script that calls the integrate_all R script

mix_config.yaml -- The configuration file for use with Snakemake on mixture files

mix_manifest_gen.R -- R script that creates a manifest file for mixture files (for use with Snakemake)

msings.sh -- Bash script used to run mSINGS in pseudobulk/aggregate expression mode (not used for manuscript results)

mtx_config.yaml -- Configure file for Snakemake to be used with MTX files

no_msi_clustering.R  -- R script used to test if removing genes used in MSIsensor-RNA baseline from samples changed clustering

patient_report_generator.R -- R script used to run Rmarkdown that creates document to plot results for each individual

plot_all_bamog.R -- R script that plots the samples for each individual (from BAM files)

plot_all_gsm.R -- R script that plots the samples for each individual (from MTX files)

pro.sh -- Bash script used to run MSIsensor-pro on pseudobulk/aggregate expression for each individual (not sued for manuscript results)

sensor2.sh -- Bash script used to run MSIsensor2 on pseudobulk/aggregate expression for each individual (not sued for manuscript results)

sensor_rna.sh -- Bash script used to run MSIsensor-RNA (also sensor_rna_shaper.R)

sensor_rna_shaper.R -- R script used to coerce single-cell data into compatible format for MSIsensor-RNA

seurat_interaction.R -- R script from inferCNV team edited to allow non-overlap in cell barcode names in integrated samples file and inferCNV object

sum_stat_catter.sh -- Bash file to quickly concatenate results from all runs (note: will overwrite any files generated from previous use without edit)

old_scripts/ -- directory containing many useful (and not useful) old versions of code or code used with the project that was deemed no longer necessary to generate the results of the manuscript (but retained for internal use)

### sensor2_results
MSIsensor2 scores for each individual in the study with a BAM file that was run at the individual and pseudobulk level (left out of final results due to lack of agreement with scATOMIC)

### sensor_rna_results
MSIsensor-rna scores for each individual in the study with a BAM file that was run at the individual and pseudobulk level (left out of final results due to lack of agreement with scATOMIC)

### summary_stats
All final results for each individual that detail the extent of heterogeneity (cluster_stats.tsv files) and more information about the ANOVA test run on the clusters (anova_results.tsv files)

Please feel free to reach out with any questions if this README has not answered your questions. 

