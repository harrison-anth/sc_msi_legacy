#!/bin/bash
module load Anaconda3
module load java

conda activate vcftools

ref=/data3/hanthony/reference_files/refdata-gex-GRCh38-2020-A/fasta/genome.fa
pon=/data3/hanthony/reference_files/protected_references/gatk4_mutect2_4136_pon.vcf.gz
gnome=/data3/hanthony/reference_files/protected_references/somatic-hg38_af-only-gnomad.hg38.vcf

while read gnorts
do
#annotate vcf files with annovar

perl ../temp/annovar/convert2annovar.pl -format vcf4 ../vcf/$gnorts.filtered.final.vcf.gz \
> ../annovar_results/$gnorts.avinput

#annotate variants
perl ../temp/annovar/annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 \
../annovar_results/$gnorts.avinput ../temp/annovar/humandb/

perl  ../temp/annovar/annotate_variation.pl -filter -dbtype exac03 -buildver hg38 \
../annovar_results/$gnorts.avinput ../temp/annovar/humandb/

#create vartrix file
conda activate seurat

vartrix -v ../vcf/$gnorts.filtered.final.vcf.gz -b ../bam/$gnorts.vcf_ready.bam -f $ref -c \
../pseudobulk_barcodes/$gnorts/just_barcodes.tsv --out-matrix ../vartrix_results/"$gnorts"_alt_matrix.mmx \
--log-level info --threads 20 \
--ref-matrix ../vartrix_results/"$gnorts"_ref_matrix.mmx \
--out-variants ../vartrix_results/"$gnorts"_var_list.tsv \
--scoring-method coverage \
--mapq 25 &

vartrix -v ../vcf/$gnorts.filtered.final.vcf.gz -b ../bam/$gnorts.vcf_ready.bam -f $ref -c \
../pseudobulk_barcodes/$gnorts/just_barcodes.tsv --out-matrix ../vartrix_results/"$gnorts"_frac_matrix.mmx \
--log-level info --threads 20 \
--scoring-method alt_frac \
--mapq 25 &

vartrix -v ../vcf/$gnorts.filtered.final.vcf.gz -b ../bam/$gnorts.vcf_ready.bam -f $ref -c \
../pseudobulk_barcodes/$gnorts/just_barcodes.tsv --out-matrix ../vartrix_results/"$gnorts"_genos_matrix.mmx \
--log-level info --threads 20 \
--scoring-method consensus \
--mapq 25 &

wait

zcat ../vcf/$gnorts.filtered.final.vcf.gz | awk '{print $1,$2}' > ../vartrix_results/$gnorts.snv_loci.txt 
sed -i 's/\s/:/g' ../vartrix_results/$gnorts.snv_loci.txt

done < ../temp/temp_vartrix_files.txt
