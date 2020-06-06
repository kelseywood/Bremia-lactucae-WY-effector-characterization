#!/bin/bash

#SBATCH -c 5  ## number of cores
#SBATCH -B 1  ## Number of machines cores are located on (1 is optimal)
#SBATCH -J PONDR ## Job name
#SBATCH -p assembly  ## partitional to submit to -- generally gc
#SBATCH --mem=10G # Memory pool required - Not currently implemented on cabernet
#SBATCH --mem-per-cpu 1000
#SBATCH --time=7-0 # day-hours:min:sec   min:sec
#SBATCH --mail-type=ALL # Tpe of email notification to recieve
#SBATCH --mail-user=mjnur@ucdavis.edu # email to which notifications will be sent. Not implemented yet.


#./runPONDR.sh NR_CLEAVED_secreted_B_lac.protein.fasta Results_secretome_cleaved_Jul29.txt

for file in RXLR_EER/*.fasta 
do
    short=$(basename "$file")
    ./runPONDR.sh $file WY_noRXLR_results_aug30/PONDR_Results_1_${short}.txt
done
