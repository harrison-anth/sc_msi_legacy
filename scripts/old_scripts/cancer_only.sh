#!/bin/bash
#SBATCH --job-name="cancer_only"
#SBATCH --output=../stdout_files/cancer_only.out
#SBATCH --error=../error_files/cancer_only.error
#SBATCH --partition="normal","highmem"
#SBATCH --mem=35G
#SBATCH --nice=555555
gnorts=$1

if [[ ! -f ../annotated_h5/"$gnorts"_cancer.rds ]]
then

module load Anaconda3
conda activate seurat

Rscript only_cancer_cells.R $gnorts

fi

