#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -t "10:00:00"

python control_population_allel_data_mbuti.py > control_mbuti.csv
python control_population_allel_data_bedouin.py > control_bedouin.csv
python control_population_allel_data_kalash.py > control_kalash.csv
python control_population_allel_data_han.py > control_han.csv
python control_population_allel_data_bougainville.py > control_bougainville.csv
python control_population_allel_data_maya.py > control_maya.csv
python control_population_allel_data_french.py > control_french.csv
