# Scripts used to generate the results of this study are located in this directory. They are described as such:

master.sh -- Used to set parameters for other scripts and select samples to be used in the analysis.

process_fastq.sh -- Run cellranger on fastq files.

b2fq.sh -- Use cellranger to convert SRA BAM file to fastq format. 

pseudobulk.sh -- Create peuduobulk data for all leiden-clusters for a sample.

msi_passer.sh -- Select MSI tools to run.




# Using snakemake
## two snakefiles exist, one of which can be used with bam files, and the other, which can be used with mtx files.
## it can be run interactively where you simply set the cores to the number of available threads:
## snakemake --use-conda --cores $(($SLURM_CPUS_ON_NODE -1))
## it can also be used with a cluster (in our case we used slurm)
## snakemake --executor cluster-generic --cluster-generic-submit-cmd "sbatch" --use-conda --latency-wait 60 --profile ../conda_envs/slurm_profile/ --verbose
## We give the cluster profile we used with this project in the ../conda_envs/slurm_profile/config.yaml file
