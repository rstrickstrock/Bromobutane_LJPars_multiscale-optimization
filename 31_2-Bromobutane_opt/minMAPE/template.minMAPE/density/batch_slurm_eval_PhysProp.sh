#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:01:00
#SBATCH --job-name=evalPhysProp

sleep 2
if [[ -f this_prediction.csv ]]; then
	break
else
	echo "MISSING: $eval_cwd/this_prediction.csv"
fi
