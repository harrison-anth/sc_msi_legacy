#!/bin/bash

if [[ -f ../mix_summary_stats/all_mix_anova_stats.tsv ]]
then
rm ../mix_summary_stats/all_mix_anova_stats.tsv
awk 'FNR>1 || NR==1' ../mix_summary_stats/*mix*anova* > ../mix_summary_stats/all_mix_anova_stats.tsv
else
awk 'FNR>1 || NR==1' ../mix_summary_stats/*mix*anova* > ../mix_summary_stats/all_mix_anova_stats.tsv
fi

if [[ -f ../mix_summary_stats/all_mix_cluster_stats.tsv ]]
then
rm ../mix_summary_stats/all_mix_cluster_stats.tsv
awk 'FNR>1 || NR==1' ../mix_summary_stats/*mix*cluster* > ../mix_summary_stats/all_mix_cluster_stats.tsv
else
awk 'FNR>1 || NR==1' ../mix_summary_stats/*mix*cluster* > ../mix_summary_stats/all_mix_cluster_stats.tsv
fi

if [[ -f ../summary_stats/all_anova_stats.tsv ]]
then
rm ../summary_stats/all_anova_stats.tsv
awk 'FNR>1 || NR==1' ../summary_stats/*anova* > ../summary_stats/all_anova_stats.tsv
else
awk 'FNR>1 || NR==1' ../summary_stats/*anova* > ../summary_stats/all_anova_stats.tsv
fi

if [[ -f ../summary_stats/all_cluster_stats.tsv ]]
then
rm ../summary_stats/all_cluster_stats.tsv 
awk 'FNR>1 || NR==1' ../summary_stats/*cluster* > ../summary_stats/all_cluster_stats.tsv
else
awk 'FNR>1 || NR==1' ../summary_stats/*cluster* > ../summary_stats/all_cluster_stats.tsv
fi

if [[ -f ../sensor_rna_results/all_cell_glob_results.txt ]]
then
rm ../sensor_rna_results/all_cell_glob_results.txt
awk 'FNR>1 || NR==1' ../sensor_rna_results/*all_cell*glob*.txt > ../sensor_rna_results/all_cell_glob_results.txt
else
awk 'FNR>1 || NR==1' ../sensor_rna_results/*all_cell*glob*.txt > ../sensor_rna_results/all_cell_glob_results.txt
fi

if [[ -f ../sensor_rna_results/cancer_cell_glob_results.txt ]]
then
rm ../sensor_rna_results/cancer_cell_glob_results.txt
awk 'FNR>1 || NR==1' ../sensor_rna_results/*cancer*glob*.txt > ../sensor_rna_results/cancer_cell_glob_results.txt
else
awk 'FNR>1 || NR==1' ../sensor_rna_results/*cancer*glob*.txt > ../sensor_rna_results/cancer_cell_glob_results.txt
fi



