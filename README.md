# Intratumoral heterogeneity in microsatellite instability status at single cell resolution

## Repository information 

This repository contains the final results and code used to generate it for our recently submitted manuscript. We will archive this github repository after the manuscript is published. 

To access the computational pipeline associated with this research, SC-MSI, please reference the following repository (https://github.com/harrison-anth/sc_msi). We chose to separate the two, as 
the legacy repository contains thousands of txt files related to the results of the manuscript, and older code not needed to run SC-MSI. 

## Description of repository structure
conda_envs -- contains all the conda environments used for the Snakemake pipeline
images -- contains UMAP and other preliminary plots for every sample
manifests -- stores all the files used to download data and lists of patients/samples
markdown_files -- Rmarkdown files used to create patient reports
mix_summary_stats -- Contains all the summary stats for each cell mixing experiment run
msings_results -- early results of samples/patients using mSINGS when we were still deciding which tool to use
pro_results -- early results of samples/patients using MSIsensor-pro when we were still deciding which tool to use
pseudobulk_barcodes -- cell barcodes for all samples and individuals
reports -- patient reports generated as part of the Snakemake pipeline
scripts -- all original code used for the study
sensor2_results -- early results of samples/patients using MSIsensor2 when we were still deciding which tool to use
sensor_rna_results -- MSIsensor-RNA results for each individual and sample (the tool we ultimately used in the study)
summary_stats -- Contains all the summary stats for each individual
