#! /bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --time=00:05:00
#SBATCH --job-name=mR2evalDens

echo -e "density\n" | gmx energy -f 1-bromobutane_prod.edr -b 100 -o 1-bromobutane_density
sleep 1
echo -e "density\n" | gmx energy -f 2-bromobutane_prod.edr -b 100 -o 2-bromobutane_density
