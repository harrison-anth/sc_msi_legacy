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
all_samples.tsv -- All samples for data from EGAD00001008555, EGAD00001008584, EGAD00001008585
sensor.yml was used to run MSIsensor MSIsensor2 and MSIsensor-pro

msings.yml was used to run mSINGS

MSIR.yml was used to run R code

mantis.yml was used to run MANTIS

vcf2maf.yml was used to convert from vcf file format to maf file format (note even with this environment, vcf2maf is notoriously finicky and will require 
a working local installation of VEP)

### Publication results
The TCGA results and scripts required to generate them are included in the tcga/ folder alongside a separate README.md
 detailing the pipeline. 

The results of all non-TCGA datasets are in the non_tcga/ folder. Each separate dataset has its own subfolder and README.md. 

### Graphs

All the code necessary to generate the graphs used in the paper are stored in the manuscript_graphs/ directory. The Rmarkdown files
used cater to my local paths and settings but can be adapated by any experienced R user. 
All graphs were generated using R version 4.1.2

### Baselines and reference files

The reference files are available from the GDC (https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)

The baselines created for each tool and used in this study are included in the baselines/ directory. 

### Test example

To ensure reproducibility of the results and quick implementation of our pipelines for other researchers who want to use these tools, we have included
a test example that allows for the user to quickly use each MSI tool on small BAM and gene count matrices that have been subset down to only the microsatellites found on 
chromosome 7. See the test/ directory for a separate README.md related to their use. 



Please feel free to reach out with any questions if this README has not answered your questions. 

