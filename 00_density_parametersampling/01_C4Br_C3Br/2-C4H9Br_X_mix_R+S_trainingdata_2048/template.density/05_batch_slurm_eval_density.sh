#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=eval_density

echo -e "density\n" | gmx energy -f 14_prod.edr -o density
