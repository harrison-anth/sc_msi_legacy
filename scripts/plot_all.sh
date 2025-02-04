#!/bin/sh
#SBATCH --job-name="plotter"
#SBATCH --output=../stdout_files/plots.out
#SBATCH --error=../error_files/plots.error
#SBATCH --partition="normal","highmem","gpu","interactive"
#SBATCH --mem=8G
#declare variables

samples=$1
gsm=$2
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/$samples)

if [[ ! -f ../images/"$gnorts".pdf ]]
then
module load Anaconda3
conda activate seurat

#Rscript plot_all.mult.R $gnorts

if [[ $gsm == "N" ]]
then

Rscript plot_all_new.R $gnorts

elif [[ $gsm == "Y" ]]
then
Rscript plot_all_gsm.R $gnorts


fi
fi

