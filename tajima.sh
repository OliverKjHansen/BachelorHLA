#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 100G
#SBATCH -c 1
#SBATCH -t "10:00:00"

python tajima_mbuti.py > tajima_mbuti.csv
python tajima_bedoiun.py > tajima_bedoiun.csv
python tajima_french.py > tajima_french.csv
python tajima_kalash.py > tajima_kalash.csv
python tajima_han.py > tajima_han.csv
python tajima_bougainville.py > tajima_bougainville.csv
python tajima_maya.py > tajima_maya.csv
python tajima_control_mbuti.py> control_tajima_mbuti.csv
python tajima_control_bedouin.py> control_tajima_bedouin.csv
python tajima_control_french.py> control_tajima_french.csv
python tajima_control_kalash.py> control_tajima_kalash.csv
python tajima_control_han.py> control_tajima_han.csv
python tajima_control_bougainville.py> control_tajima_bougainville.csv
python tajima_control_maya.py> control_tajima_maya.csv
