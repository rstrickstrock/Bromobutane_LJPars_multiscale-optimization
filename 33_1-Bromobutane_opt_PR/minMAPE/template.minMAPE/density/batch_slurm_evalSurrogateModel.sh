#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalDensSurrModel

python loadEvalSurrogModel.py

