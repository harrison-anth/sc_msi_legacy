#!/bin/sh
#SBATCH --job-name="anno"
#SBATCH --output=../stdout_files/anno.out
#SBATCH --error=../error_files/anno.error
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=8G
#declare variables

samples=$1
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/$samples)
gsm=$2

if [[ ! -f ../annotated_h5/"$gnorts".rds ]]
then

module load Anaconda3
conda activate seurat

if [[ $gsm == "Y" ]]
then

Rscript annotate_gsmog.R $gnorts

elif [[ $gsm = "N" ]]
Rscript annotate_bamog.R $gnorts

fi
fi
