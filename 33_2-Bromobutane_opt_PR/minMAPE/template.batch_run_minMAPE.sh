#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=24:00:00
#SBATCH --job-name=2-mM-QQ

./run_minMAPE-QQ.sh
