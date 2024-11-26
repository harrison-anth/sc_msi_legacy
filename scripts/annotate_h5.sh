#!/bin/sh
#SBATCH --job-name="anno"
#SBATCH --output=../stdout_files/anno.out
#SBATCH --error=../error_files/anno.error
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=8G
#declare variables

samples=$1
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/$samples)
flavor=$2

if [[ "$flavor" == "depth" ]]
then
if [[ ! -f ../annotated_h5/"$gnorts".rds ]]
then
module load Anaconda3
conda activate seurat
Rscript annotate_h5.fast.new_mut.R $gnorts
fi
elif [[ "$flavor" == "alt" ]]
then
if [[ ! -f ../annotated_h5/"$gnorts"_alt_count.rds ]]
then
module load Anaconda3
conda activate seurat
Rscript annotate_h5.fast.1.R $gnorts
fi
elif [[ "$flavor" == "af" ]]
then
if [[ ! -f ../annotated_h5/"$gnorts"_af.rds ]]
then
module load Anaconda3
conda activate seurat
Rscript annotate_h5.fast.2.R $gnorts
fi
fi




