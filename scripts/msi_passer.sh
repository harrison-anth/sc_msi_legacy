#!/bin/bash
#SBATCH --output=../stdout_files/msi_passer.out
#SBATCH --error=../error_files/msi_passer.error
#SBATCH --partition="normal","highmem","gpu"

gnorts=$1
gsm=$2

#which tools
pro=Y
msings=Y
sensor2=Y
premsim=Y


#make vcfs?
vcf=N


count=$(ls ../bam/"$gnorts"/*.bam | wc -l )

if [[ "$sensor2" == "Y" -a "$gsm" == "N"  ]]
then
sbatch --array 1-$count sensor2.sh $gnorts
fi

if [[ $msings == "Y" -a "$gsm" == "N" ]]
then
sbatch --array 1-$count msings.sh $gnorts
fi

if [[ $pro == "Y" -a "$gsm" == "N" ]]
then
sbatch --array 1-$count pro.sh $gnorts
fi

if [[ $vcf == "Y" -a "$gsm" == "N" ]]
then 
sbatch --array 1-$count make_vcf.sh $gnorts
fi

if [[ $premsim == "Y" ]]
then
sbatch sensor_rna.sh $gnorts
fi


