#!/bin/sh
#SBATCH --job-name="sc_pro"
#SBATCH --output=../stdout_files/pro.out
#SBATCH --error=../error_files/pro.error
#SBATCH --partition="normal","highmem","gpu"

gnorts=$1
num=$((SLURM_ARRAY_TASK_ID-1))
baseline=/data4/hanthony/tcga_msi_tools/baselines/pro_rna/scan_tcga_ms.bed_baseline
tumor=../bam/$gnorts/$num.bam

if [[ ! -f ../pro_results/"$gnorts"_cluster_"$num" ]]
then
module load Anaconda3
conda activate pro
msisensor-pro pro -d $baseline -t $tumor -o ../pro_results/"$gnorts"_cluster_"$num"
rm ../pro_results/"$gnorts"_cluster_"$num"_unstable
rm ../pro_results/"$gnorts"_cluster_"$num"_all
rm ../pro_results/"$gnorts"_cluster_"$num"_dis

fi

