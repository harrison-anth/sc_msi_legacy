#!/bin/sh
#SBATCH --job-name="vcf_sin"
#SBATCH --output=../stdout_files/vcf.out
#SBATCH --error=../error_files/vcf.error
#SBATCH --partition="normal","highmem","gpu"
#SBATCH --mem=25G
#SBATCH --cpus-per-task=16
#SBATCH --nice=55555

module load Anaconda3
module load java

ref=/data3/hanthony/reference_files/refdata-gex-GRCh38-2020-A/fasta/genome.fa
pon=/data3/hanthony/reference_files/protected_references/gatk4_mutect2_4136_pon.vcf.gz
gnome=/data3/hanthony/reference_files/protected_references/somatic-hg38_af-only-gnomad.hg38.vcf



gnorts=$1
num=$((SLURM_ARRAY_TASK_ID-1))
tum_id=$gnorts

if [[ ! -f ../bam/$gnorts/$num.vcf_ready.bam ]]
then
tumor=../bam/$gnorts/$num.bam
module load Anaconda3
module load java
conda activate vcftools

gatk SplitNCigarReads \
-R $ref \
-I $tumor \
-O ../vcf/$gnorts.$num.vcf_ready.bam

mv ../vcf/$gnorts."$num".vcf_ready.bai ../vcf/$gnorts."$num".vcf_ready.bam.bai



fi

tumor=../vcf/$gnorts.$num.vcf_ready.bam

if [[ ! -f ../vcf/$tum_id.$num.vcf.gz ]]
then
gatk Mutect2 \
-R $ref \
-I $tumor \
--panel-of-normals $pon \
--native-pair-hmm-threads $SLURM_CPUS_PER_TASK \
-O ../vcf/$tum_id.$num.vcf.gz
fi

if [[ ! -f ../vcf/$tum_id.$num.filtered.vcf.gz ]]
then
conda activate vcftools

gatk FilterMutectCalls \
-V ../vcf/$tum_id.$num.vcf.gz \
-R $ref \
-O ../vcf/$tum_id.$num.filtered.vcf.gz
fi

#annotate vcf files with annovar

#convert vcf to annovar format
if [[ ! -f ../annovar_results/$gnorts.$num.avinput ]] 
then
conda activate vcftools
perl ../temp/annovar/convert2annovar.pl -format vcf4 ../vcf/$gnorts.$num.filtered.vcf.gz \
> ../annovar_results/$gnorts.$num.avinput
fi

#annotate variants
if [[ ! -f ../annovar_results/$gnorts.$num.avinput.exonic_variant_function ]]
then
conda activate vcftools
perl ../temp/annovar/annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 \
../annovar_results/$gnorts.$num.avinput ../temp/annovar/humandb/
fi

#get cytoband info
if [[ ! -f ../annovar_results/$gnorts.$num.avinput.hg38_cytoBand ]]
then
conda activate vcftools
perl ../temp/annovar/annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg38 \
../annovar_results/$gnorts.$num.avinput ../temp/annovar/humandb/
fi

#create filtered variant info
if [[ ! -f ../annovar_results/$gnorts.$num.avinput.hg38_exac03_filtered ]]
then 
conda activate vcftools
perl  ../temp/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 \
../annovar_results/$gnorts.$num.avinput ../temp/annovar/humandb/
fi

#create vartrix file

if [[ ! -f ../vartrix_results/$gnorts.$num.mmx ]]
then
conda activate seurat

vartrix -v ../vcf/$tum_id.$num.filtered.vcf.gz -b ../bam/$gnorts/$num.vcf_ready.bam -f $ref -c \
../pseudobulk_barcodes/$gnorts/"$gnorts"_cluster_"$num".tsv -o ../vartrix_results/$gnorts.$num.mmx
fi

if [[ ! -f ../vartrix_results/$gnorts.$num.snv_loci.txt ]]
then
zcat ../vcf/$gnorts.$num.filtered.vcf.gz | awk '{print $1,$2}' > ../vartrix_results/$gnorts.$num.snv_loci.txt 
sed -i 's/\s/:/g' ../vartrix_results/$gnorts.$num.snv_loci.txt
fi

if [[ ! -f ../vartrix_results/$gnorts.$num.rds ]]
then
module load Anaconda3
conda activate seurat
Rscript gedder.R $gnorts $num
fi

if [[ ! -f ../maf/"$gnorts"/$gnorts.$num.maf ]]
then
if [[ ! -d ../maf/"$gnorts"/ ]]
then
mkdir ../maf/"$gnorts"/
fi

module load Anaconda3
conda activate vcf22maf

vcf=../vcf/$gnorts.$num.filtered.vcf.gz
gunzip $vcf
vcf=../vcf/$gnorts.$num.filtered.vcf

vcf2maf.pl --input $vcf --output-maf ../maf/"$gnorts"/$gnorts.$num.maf \
--ref-fasta $ref --tumor-id $gnorts.$SLURM_ARRAY_TASK_ID --vep-path /data3/hanthony/.conda_envs/vcf22maf/bin \
--vep-data /data4/hanthony/vep_files/ --ncbi-build GRCh38 --vep-forks 1

gzip $vcf

fi
