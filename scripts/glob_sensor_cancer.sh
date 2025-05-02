#!/bin/sh

module load Anaconda3

while read patient_id
do


if [[ ! -f ../sensor_rna_results/"$patient_id"_cancer_glob.txt ]]
then


conda activate seurat


#create global pseudobulk expression for sample

Rscript glob_shaper_cancer.R $patient_id


conda deactivate 

data=../temp/"$patient_id"_glob_cancer_training_model.csv

new_model=../temp/"$patient_id"_cancer_glob.model

msisensor-rna train -i $data -m $new_model -t PanCancer

#global pseudobulk expression results
expr=../temp/"$patient_id"_cancer_glob.csv


msisensor-rna detection -i $expr -o ../sensor_rna_results/"$patient_id"_cancer_glob.txt -m $new_model -d True

fi

done < ../manifests/glob_patients.txt
