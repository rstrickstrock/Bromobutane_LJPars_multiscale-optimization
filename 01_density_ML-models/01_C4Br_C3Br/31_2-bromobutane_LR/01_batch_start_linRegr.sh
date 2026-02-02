#!/bin/bash

#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=10:00:00
#SBATCH --job-name=linRegr

python 01_linearRegression.py

