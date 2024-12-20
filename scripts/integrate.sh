#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --error=../error_files/integrate.error
#SBATCH --output=../stdout_files/integrate.out
#SBATCH --partition="normal","highmem","gpu","interactive"
#SBATCH --mem=16G

#define variables
gsm=$1




if [[ $gsm == "Y" ]]
then
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/gsm_patients.txt)
else
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/patient_ids.csv)
fi


if [[ ! -f ../integrated_samples/"$gnorts".rds ]]
then

module load Anaconda3
conda activate atomic

Rscript integrate.R $gnorts $gsm

fi


if [[ -f ../integrated_samples/"$gnorts".rds ]]
then

module load Anaconda3
conda activate seurat

Rscript plot_integrate.R $gnorts $gsm

fi


#for pseudobulk sensor-rna analysis

sbatch sensor_rna.sh $gnorts

sbatch sensor_rna.sh "$gnorts"_cancer

#generate report; will only work if above sensor-rna scripts were successful 

if [[ ! -f ../reports/$gnorts.html ]]
then
module load Anaconda3
conda activate atomic

Rscript patient_report_generator.R $gnorts
 
fi





