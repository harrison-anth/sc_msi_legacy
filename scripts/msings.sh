#!/bin/sh
#SBATCH --job-name="sc_msings"
#SBATCH --output=../stdout_files/msings.out
#SBATCH --error=../error_files/msings.error
<<<<<<< HEAD
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=20G
=======
#SBATCH --partition="normal","highmem"
#SBATCH --mem=20G

module load Anaconda3

conda activate msings
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2


gnorts=$1

num=$((SLURM_ARRAY_TASK_ID-1))
MSI_BASELINE=/data4/hanthony/tcga_msi_tools/baselines/msings_50_filtered.baseline
tumor=../bam/$gnorts/$num.bam
BEDFILE=/data4/hanthony/tcga_msi_tools/bed_files/ms_sites_premium+.bed
REF_GENOME=/data3/hanthony/reference_files/hg38.fa
msings=/home/hanthony/bin/programs/msings/msi


<<<<<<< HEAD
=======
if [[ ! -f $tumor.bai ]]
then
samtools index $tumor
fi

>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
if [[ ! -f ../msings_results/"$gnorts"_cluster_"$num".MSI_Analysis.txt ]]
then
module load Anaconda3
conda activate msings

mkdir -p ../temp/"$gnorts"_"$num"
outdir=../temp/"$gnorts"_"$num"

<<<<<<< HEAD
=======
mkdir -p ../temp/"$gnorts"_"$num"
outdir=../temp/"$gnorts"_"$num"

>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
samtools mpileup -f $REF_GENOME -d 100000 -A -E -l $BEDFILE $tumor | awk '{if($4 >= 6) print $0}' > $outdir/$gnorts.mpileup

$msings analyzer $outdir/$gnorts.mpileup $BEDFILE -o $outdir/$gnorts.msi.txt

$msings count_msi_samples $MSI_BASELINE $outdir -m $multiplier -t $msi_min_threshold $msi_max_threshold \
-o ../msings_results/"$gnorts"_cluster_"$num".MSI_Analysis.txt

if [[ -f ../msings_results/"$gnorts"_cluster_"$num".MSI_Analysis.txt ]]
then
#cleanup
rm -r $outdir
fi

fi
