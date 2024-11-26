#!/bin/sh
#SBATCH --job-name="bash"
#SBATCH --output=../stdout_files/anno.out
#SBATCH --error=../error_files/anno.error
#SBATCH --partition="interactive"
#SBATCH --nodes=1
#SBATCH --nice=55555
#declare variables

samples=$1
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/$samples)

if [[ ! -f ../annotated_h5/"$gnorts".rds ]]
then
module load Anaconda3
conda activate seurat

Rscript annotate_h5.fast.debug.R $gnorts

fi


