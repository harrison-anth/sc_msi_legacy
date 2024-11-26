#!/bin/bash
#SBATCH --job-name=numbat
#SBATCH --error=../error_files/numbat.error
#SBATCH --output=../stdout_files/numbat.out
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=30G

#define variables
gnorts=$1
numbat=/data3/hanthony/.conda_envs/seurat/lib/R/library/numbat/bin/pileup_and_phase.R

SLURM_CPUS_PER_TASK=1

if [[ ! -f ../numbat/"$gnorts"/"$gnorts"_allele_counts.tsv.gz ]]
then

if [[ ! -d ../numbat/"$gnorts"/ ]]
then
mkdir ../numbat/"$gnorts"/
fi

module load Anaconda3
conda activate seurat

barcodes=../pseudobulk_barcodes/"$gnorts"/just_barcodes.tsv
eagle=/data4/hanthony/single_msi/temp/Eagle_v2.4.1/eagle
ggmap=/data4/hanthony/single_msi/temp/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
vcf=/data4/hanthony/single_msi/temp/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
panel=/data4/hanthony/single_msi/temp/1000G_hg38


Rscript $numbat --samples $gnorts --bams ../cell_ranger_output/"$gnorts"/outs/possorted_genome_bam.bam \
--barcodes $barcodes --gmap $ggmap --eagle $eagle --snpvcf $vcf --paneldir $panel \
--ncores $SLURM_CPUS_PER_TASK --outdir ../numbat/"$gnorts"/

fi


if [[ ! -f ../numbat/"$gnorts"/bulk_clones_final.tsv.gz ]]
then
module load Anaconda3
conda activate seurat

Rscript numbat.R $gnorts $SLURM_CPUS_PER_TASK

fi


