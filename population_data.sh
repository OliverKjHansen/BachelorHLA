#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -t "10:00:00"

python population_allel_data_mbuti.py > mbuti.csv
python population_allel_data_bedouin.py > bedouin.csv
python population_allel_data_kalash.py > kalash.csv
python population_allel_data_han.py > han.csv
python population_allel_data_bougainville.py > bougainville.csv
python population_allel_data_maya.py > maya.csv
python population_allel_data_french.py > french.csv