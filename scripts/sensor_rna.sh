#!/bin/sh
#SBATCH --job-name="premsim"
#SBATCH --output=../stdout_files/premsim.out
#SBATCH --error=../error_files/premsim.error
#SBATCH --partition="normal","highmem","gpu"

gnorts=$1

module load Anaconda3

#if [[ ! -f ../premsim/$gnorts.rds ]]
#then

#conda activate seurat
#
#Rscript sensor_rna_shaper.R $gnorts
#
#fi

if [[ ! -f ../sensor_rna_results/$gnorts.txt ]]
then

conda deactivate

data=../temp/"$gnorts"_training_model.csv

msisensor-rna train -i ../temp/"$gnorts"_training_model.csv -m ../temp/"$gnorts".model -t PanCancer

msisensor-rna detection -i ../temp/$gnorts.csv -o ../sensor_rna_results/$gnorts.txt -m ../temp/$gnorts.model -d True

fi
