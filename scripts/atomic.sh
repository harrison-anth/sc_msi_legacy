#!/bin/bash
#SBATCH --job-name=atomic
#SBATCH --error=../error_files/atomic.error
#SBATCH --output=../stdout_files/atomic.out
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --nodes=1

#define variables
gnorts=$1



if [[ ! -f ../atomic/"$gnorts".rds ]]
then

module load Anaconda3
conda activate atomic

Rscript atomic.R $gnorts 

fi

#sbatch cancer_only.sh $gnorts
