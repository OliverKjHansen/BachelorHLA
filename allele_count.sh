#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -t "10:00:00"
 
python allele_count.py > allel_count.csv
python allele_count_control.py > allel_count_control.csv