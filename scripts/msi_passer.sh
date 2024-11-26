#!/bin/bash
#SBATCH --output=../stdout_files/msi_passer.out
#SBATCH --error=../error_files/msi_passer.error
#SBATCH --partition="normal","highmem","gpu"

gnorts=$1

#which tools
<<<<<<< HEAD
pro=N
msings=N
sensor2=N
premsim=Y


#make vcfs?
vcf=N

=======
pro=Y
msings=Y
sensor2=Y
>>>>>>> 8a5c77b783c5b63b0d7691a18aacc16be0338eef

count=$(ls ../bam/"$gnorts"/*.bam | wc -l )

if [[ $sensor2 == "Y" ]]
then
sbatch --array 1-$count sensor2.sh $gnorts
fi

if [[ $msings == "Y" ]]
then
sbatch --array 1-$count msings.sh $gnorts
fi

if [[ $pro == "Y" ]]
then
sbatch --array 1-$count pro.sh $gnorts
fi

if [[ $vcf == "Y" ]]
then 
sbatch --array 1-$count make_vcf.sh $gnorts
fi

if [[ $premsim == "Y" ]]
then
sbatch sensor_rna.sh $gnorts
fi


