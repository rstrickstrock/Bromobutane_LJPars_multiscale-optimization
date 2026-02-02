#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=48:00:00
#SBATCH --job-name=GPR-6

python 01_gaussianProcessRegression-6.py
