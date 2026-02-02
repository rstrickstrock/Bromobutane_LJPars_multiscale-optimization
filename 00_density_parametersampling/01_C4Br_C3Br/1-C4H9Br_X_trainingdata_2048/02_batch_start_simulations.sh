#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=72:00:00
#SBATCH --job-name=1BromoDens

python 02_start_simulations.py
