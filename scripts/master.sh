#!/bin/bash

#set to Y to annotate exisiting
first_run=N

<<<<<<< HEAD
annotate_mode=N
=======
annotate_mode=Y
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef

plot_mode=N

integrate_mode=Y

integrate_cancer=N

gsm=N
	
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
<<<<<<< HEAD
then 
sbatch --array 1-$count process_fastq.sh $samples $gsm
=======
then sbatch --array 10-$count%10 process_fastq.sh $samples
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef
fi

#if [[ $annotate_mode == "Y" ]]
#then
#flavor=af
#sbatch --array 1-$count annotate_h5.sh $samples $flavor
#fi

if [[ $annotate_mode == "Y" ]]
then
flavor=depth
sbatch --array 1-$count annotate_h5.sh $samples $flavor
fi

#if [[ $annotate_mode == "Y" ]]
##then
#flavor=alt
#sbatch --array 1-$count annotate_h5.sh $samples $flavor
#fi

if [[ $plot_mode == "Y" ]]
then

sbatch --array 1-$count plot_all.sh $samples
fi

if [[ $integrate_mode == "Y" ]]
then
sbatch --array 1-18 integrate.sh $gsm
fi



if [[ $integrate_cancer == "Y" ]]
then
sbatch --array 1-18 integrate_cancer.sh
fi
