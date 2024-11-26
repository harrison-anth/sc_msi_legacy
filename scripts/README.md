# Scripts used to generate the results of this study are located in this directory. They are described as such:

master.sh -- Used to set parameters for other scripts and select samples to be used in the analysis.

process_fastq.sh -- Run cellranger on fastq files.

b2fq.sh -- Use cellranger to convert SRA BAM file to fastq format. 

pseudobulk.sh -- Create peuduobulk data for all leiden-clusters for a sample.

msi_passer.sh -- Select MSI tools to run.
