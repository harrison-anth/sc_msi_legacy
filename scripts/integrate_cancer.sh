#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --error=../error_files/integrate.error
#SBATCH --output=../stdout_files/integrate.out
#SBATCH --partition="normal","highmem","gpu","interactive"
#SBATCH --mem=8GB

#define variables


gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/patient_ids.csv)

if [[ ! -f ../integrated_samples/"$gnorts".rds ]]
then

module load Anaconda3
conda activate atomic

Rscript integrate_cancer.R $gnorts 

fi


