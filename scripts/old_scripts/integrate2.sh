#!/bin/bash
#SBATCH --job-name=integrate2
#SBATCH --error=../error_files/integrate2.error
#SBATCH --output=../stdout_files/integrate2.out
#SBATCH --partition="normal","highmem","gpu","interactive"
#SBATCH --mem=15G

#define variables
gnorts=$1
gsm=$2

if [[ ! -f ../integrated_samples/"$gnorts"2.rds ]]
then

module load Anaconda3
conda activate atomic

Rscript integrate2.R $gnorts $gsm

fi


if [[ -f ../integrated_samples/"$gnorts".rds ]]
then

module load Anaconda3
conda activate seurat

Rscript plot_integrate2.R $gnorts $gsm

fi
