#!/bin/bash 
#SBATCH -c 2  ## number of cores
#SBATCH -B 1  ## Number of machines cores are located on (1 is optimal)
#SBATCH -J ALL ## Job name
#SBATCH -p production ## partitional to submit to -- generally gc
#SBATCH --mem=8G # Memory pool required - Not currently implemented on cabernet
#SBATCH --mem-per-cpu 1000
#SBATCH --time=7-0 # day-hours:min:sec   min:sec
#SBATCH --mail-type=END # Tpe of email notification to recieve
#SBATCH --mail-user=mjnur@ucdavis.edu # email to which notifications will be sent. Not implemented yet.

for file in $(ls /share/rwmwork/mjnur/effectR/motif_redo_201810/genomes/NR_cleaved_secretomes/*.fasta | sort)
do
	absPath=$(readlink -e $file)
	echo "searching $absPath"
	Rscript generate_motif_table.R $absPath

done

