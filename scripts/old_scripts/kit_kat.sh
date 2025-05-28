#!/bin/bash
#SBATCH --partition="normal","highmem"
#SBATCH --job-name=copykat
#SBATCH --error=../error_files/copykat.error
#SBATCH --output=../stdout_files/copykat.out
#SBATCH --partition="normal","highmem"
#SBATCH --mem=30G
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1


gnorts=$1
gsm=$2

if [[ ! -f ../copy_kat_results/"$gnorts"_copykat_prediction.txt ]]
then

module load Anaconda3
conda activate MSIR

Rscript kit_kat.R $gnorts $SLURM_CPUS_PER_TASK $gsm

mv "$gnorts"_copykat* ../copy_kat_results/

fi
