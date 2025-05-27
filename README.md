# Performance assessment of computational tools to detect microsatellite instability
## Article DOI: https://doi.org/10.1093/bib/bbae390

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie

## Repository information

This repository contains all the code needed to reproduce the benchmarking results of our recent manuscript.
We compared the performance of several MSI tools on different next-generation sequencing datasets.
Examples of how to run each MSI tool are provided in the test/ directory. The structure of the GitHub repository is described
below.

### Conda environments
The conda_envs folder contains all conda environments used for data handling and all MSI tools.

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

