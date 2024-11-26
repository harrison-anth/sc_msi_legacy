#!/bin/sh
#SBATCH --job-name="handle FASTQ"
#SBATCH --output=../stdout_files/fastq%.out
#SBATCH --error=../error_files/fastq%.error
#SBATCH --partition="normal","highmem"
<<<<<<< HEAD
=======
#SBATCH --exclusive
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
#SBATCH --mem=30G

#declare variables

samples=$1
<<<<<<< HEAD
gsm=$2

=======
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
gnorts=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../manifests/$samples)
cellranger=/home/hanthony/bin/programs/cellranger-7.2.0/cellranger
transcriptome=/data3/hanthony/reference_files/refdata-gex-GRCh38-2020-A

module load Anaconda3
<<<<<<< HEAD

#make sure cellranger ran successfully

if [[ $gsm == "N" ]]
then

sanity=Y
if [[ $sanity == "Y" ]]
then
if [[ ! -d ../cell_ranger_output/"$gnorts"/outs ]]
then
rm -r ../cell_ranger_output/"$gnorts"/
fi
fi

=======
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef


if [[ ! -d ../cell_ranger_output/"$gnorts" ]]
then
mkdir ../cell_ranger_output/"$gnorts"
$cellranger count --id=$gnorts --transcriptome=$transcriptome --fastqs=../fastq/ --sample=$gnorts \
--output-dir ../cell_ranger_output/$gnorts
fi

if [[ ! -f ../sensor2_results/"$gnorts"_msi_status ]]
then
conda activate pro

tumor=../cell_ranger_output/"$gnorts"/outs/possorted_genome_bam.bam
models=/data4/hanthony/tcga_msi_tools/baselines/hg38_models

msisensor2 msi -d /data3/hanthony/reference_files/scan_tcga_ms.bed -t $tumor \
-M $models -c 5 -o ../sensor2_results/"$gnorts"_msi_status
fi

if [[ ! -f ../pro_results/"$gnorts"_msi_status ]]
then
conda activate pro
<<<<<<< HEAD
baseline=/data4/hanthony/tcga_msi_tools/baselines/pro_rna/scan_tcga_ms.bed_baseline
=======
baseline=../temp/pro_test_baseline.txt
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
tumor=../cell_ranger_output/"$gnorts"/outs/possorted_genome_bam.bam

msisensor-pro pro -d $baseline -t $tumor -o ../pro_results/"$gnorts"_msi_status

fi



if [[ ! -f ../msings_results/"$gnorts"_msi_status.MSI_Analysis.txt ]]
then
conda activate msings
<<<<<<< HEAD

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2
MSI_BASELINE=/data4/hanthony/tcga_msi_tools/baselines/msings_50_filtered.baseline
tumor=../cell_ranger_output/"$gnorts"/outs/possorted_genome_bam.bam
BEDFILE=/data4/hanthony/tcga_msi_tools/bed_files/ms_sites_premium+.bed
REF_GENOME=/data3/hanthony/reference_files/hg38.fa
msings=/home/hanthony/bin/programs/msings/msi


mkdir -p ../temp/"$gnorts"_dir
outdir=../temp/"$gnorts"_dir
samtools mpileup -f $REF_GENOME -d 100000 -A -E -l $BEDFILE $tumor | awk '{if($4 >= 6) print $0}' > $outdir/$gnorts.mpileup

$msings analyzer $outdir/$gnorts.mpileup $BEDFILE -o $outdir/$gnorts.msi.txt

$msings count_msi_samples $MSI_BASELINE $outdir -m $multiplier -t $msi_min_threshold $msi_max_threshold \
-o ../msings_results/"$gnorts"_msi_status.MSI_Analysis.txt

if [[ -f ../msings_results/"$gnorts"_msi_status.MSI_Analysis.txt ]]
then
#cleanup
rm -r $outdir
fi

fi


fi

sbatch pseudobulk.sh $gnorts $gsm
=======

#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2
MSI_BASELINE=/data4/hanthony/tcga_msi_tools/baselines/msings_rna_filtered.baseline
tumor=../cell_ranger_output/"$gnorts"/outs/possorted_genome_bam.bam
BEDFILE=/data4/hanthony/tcga_msi_tools/bed_files/ms_sites_premium+.bed
REF_GENOME=/data3/hanthony/reference_files/hg38.fa
msings=/home/hanthony/bin/programs/msings/msi


mkdir -p ../temp/"$gnorts"_dir
outdir=../temp/"$gnorts"_dir
samtools mpileup -f $REF_GENOME -d 100000 -A -E -l $BEDFILE $tumor | awk '{if($4 >= 6) print $0}' > $outdir/$gnorts.mpileup

$msings analyzer $outdir/$gnorts.mpileup $BEDFILE -o $outdir/$gnorts.msi.txt

$msings count_msi_samples $MSI_BASELINE $outdir -m $multiplier -t $msi_min_threshold $msi_max_threshold \
-o ../msings_results/"$gnorts"_msi_status.MSI_Analysis.txt

if [[ -f ../msings_results/"$gnorts"_msi_status.MSI_Analysis.txt ]]
then
#cleanup
rm -r $outdir
fi

fi




sbatch pseudobulk.sh $gnorts
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef




