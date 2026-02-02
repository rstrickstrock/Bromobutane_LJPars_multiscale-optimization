#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 32G
#SBATCH --time=15:00:00
#SBATCH --job-name=NN-14

python 01_NNTraining-14.py
