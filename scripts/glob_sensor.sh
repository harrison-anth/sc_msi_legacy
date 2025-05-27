#!/bin/sh


while read patient_id
do


#create global pseudobulk expression for sample

if [[ ! -f ../sensor_rna_results/"$patient_id"_all_cell_glob.txt ]]
then

conda activate seurat
Rscript glob_shaper.R $patient_id
fi


if [[ ! -f ../sensor_rna_results/"$patient_id"_all_cell_glob.txt ]]
then


conda deactivate 
data=../temp/"$patient_id"_glob_training_model.csv
new_model=../temp/"$patient_id"_glob.model
msisensor-rna train -i $data -m $new_model -t PanCancer
#global pseudobulk expression results
expr=../temp/"$patient_id"_glob.csv


msisensor-rna detection -i $expr -o ../sensor_rna_results/"$patient_id"_all_cell_glob.txt -m $new_model -d True


fi

done < ../manifests/all_patients.txt

