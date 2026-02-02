#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=48:00:00
#SBATCH --job-name=GPR-9

python 01_gaussianProcessRegression-9.py
