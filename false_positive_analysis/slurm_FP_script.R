#!/bin/bash

#SBATCH -c 2  ## number of cores
#SBATCH -B 1  ## Number of machines cores are located on (1 is optimal)
#SBATCH -J FP ## Job name
#SBATCH -p production  ## partitional to submit to -- generally gc
#SBATCH --mem=3G # Memory pool required - Not currently implemented on cabernet
#SBATCH --mem-per-cpu 1000
#SBATCH --time=7-0 # day-hours:min:sec   min:sec
#SBATCH --mail-type=ALL # Tpe of email notification to recieve
#SBATCH --mail-user=mjnur@ucdavis.edu # email to which notifications will be sent. Not implemented yet.

Rscript FP_script.R
