#!/bin/sh
#SBATCH --job-name="b2fq"
#SBATCH --output=../stdout_files/b2fq.out
#SBATCH --error=../error_files/b2fq.error
#SBATCH --partition="normal","highmem"
#SBATCH --array=1-6

cellranger=/home/hanthony/bin/programs/cellranger-7.2.0/cellranger
ref=/data4/hanthony/sc_rna_msi/references/refdata-gex-GRCh38-2020-A
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/6bam.txt)



$cellranger bamtofastq ../bam/$gnorts.sorted.bam ../fastq/
