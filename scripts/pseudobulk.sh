#!/bin/bash
#SBATCH --job-name="pseudobulk"
#SBATCH --output=../stdout_files/pseudobulk.out
#SBATCH --error=../error_files/pseudobulk.error
#SBATCH --partition="normal","highmem"
#SBATCH --mem=35G

gnorts=$1
gsm=$2

#run numbat?
numbat=Y
#run copykat?
copykat=Y
#run atomic?
atomic=Y
#run msi tools?
msi=Y


if [[ ! -d ../pseudobulk_barcodes/$gnorts ]]
then
module load Anaconda3
conda activate seurat

mkdir ../pseudobulk_barcodes/$gnorts

Rscript barcode_generator.R $gnorts $gsm

fi

if [[ $gsm == "N" ]]
then

if [[ ! -d ../bam/$gnorts/ ]]
then
module load Anaconda3
conda activate rusty

mkdir ../bam/$gnorts/

sinto filterbarcodes -b ../cell_ranger_output/$gnorts/outs/possorted_genome_bam.bam \
-c ../pseudobulk_barcodes/"$gnorts"/"$gnorts"_all_cell_barcodes.tsv \
--barcodetag "CB" -p 20 --outdir ../bam/$gnorts/

i=1 

for bam in $(ls -v ../bam/$gnorts/*.bam)
do

samtools view -H $bam  | sed "s/SM:[^\t]*/SM:$gnorts.$i/g" | samtools reheader - $bam > $bam.bam
mv $bam.bam $bam

samtools index $bam



i=$(( $i+1 ))

done

fi

fi

if [[ ! -f ../pseudobulk_barcodes/$gnorts/just_barcodes.tsv ]]
then
cut -f 1 ../pseudobulk_barcodes/$gnorts/"$gnorts"_all_cell_barcodes.tsv > ../pseudobulk_barcodes/$gnorts/just_barcodes.tsv
fi



if [[ $msi == "Y" ]]
then
sbatch msi_passer.sh $gnorts $gsm
fi
if [[ $copykat == "Y" ]]
then
sbatch kit_kat.sh $gnorts $gsm
fi
if [[ $numbat == "Y" ]]
then
sbatch numbat.sh $gnorts
fi
if [[ $atomic == "Y" ]]
then
sbatch atomic.sh $gnorts
fi

