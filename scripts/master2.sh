#!/bin/bash

#set to Y to annotate exisiting
first_run=N

#set sample batch
gsm=Y


annotate_mode=N
	
plot_mode=N


integrate_mode=Y

integrate_cancer=N




if [[ $gsm == "Y" ]]
then
#geo expression samples with treatment metadata
samples=all_gsm_samples.txt
gsm=Y
else
#SRA samples with MSI-H clinical metadata
samples=all_samples.tsv
gsm=N

fi 

count=$(echo $(( $(wc -l < ../manifests/$samples) )))

#flavor=af
#flavor=depth
#flavor=alt


if [[ $first_run == "Y" ]]
then 
sbatch --array 1-$count process_fastq.sh $samples $gsm
fi

if [[ $annotate_mode == "Y" ]]
then

sbatch --array 1-$count annotate_h5_new.sh $samples $gsm

fi


if [[ $plot_mode == "Y" ]]
then

sbatch --array 1-$count plot_all.sh $samples $gsm

fi

if [[ $integrate_mode == "Y" ]]
then
sbatch --array 1-18 integrate.sh $gsm
fi

