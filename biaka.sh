#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -t "10:00:00"

python population_allel_data_biaka.py > biaka.csv