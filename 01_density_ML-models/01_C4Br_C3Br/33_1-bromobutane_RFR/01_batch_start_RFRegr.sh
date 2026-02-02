#!/bin/bash

#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=24:00:00
#SBATCH --job-name=RFG

python 01_randomForestRegression.py

