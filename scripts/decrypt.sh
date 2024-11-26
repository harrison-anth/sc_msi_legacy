#!/bin/bash

module load Anaconda3
conda activate aws

C4GH_PASSPHRASE=Jet575

while read line
do
crypt4gh decrypt --sk /home/hanthony/.ssh/ega < ../encrypted/prod_ega-box-2056_"$line".*.c4gh > ../decrypted/$line

done < ../manifests/just_file_names.txt
