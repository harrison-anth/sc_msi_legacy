#!/bin/sh
#SBATCH --job-name="sensor2"
#SBATCH --output=../stdout_files/sensor2.out
#SBATCH --error=../error_files/sensor2.error
#SBATCH --partition="normal","highmem","gpu"

gnorts=$1

num=$((SLURM_ARRAY_TASK_ID-1))

models=/data4/hanthony/tcga_msi_tools/baselines/hg38_models
tumor=../bam/$gnorts/$num.bam

if [[ ! -f ../sensor2_results/"$gnorts"_cluster_"$num" ]]
then
module load Anaconda3
conda activate pro
msisensor2 msi -d /data3/hanthony/reference_files/scan_tcga_ms.bed -t $tumor -M $models -c 5 -o ../sensor2_results/"$gnorts"_cluster_"$num"

fi
