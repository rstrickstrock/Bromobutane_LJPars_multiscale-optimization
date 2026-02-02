#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalRCE


echo "Total-Energy" | gmx energy -f output/21_nvt.edr -o output/21_energy >>output/21_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/22_nvt.edr -o output/22_energy >>output/22_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/23_nvt.edr -o output/23_energy >>output/23_totalEnergy.txt
sleep 2

python gatherResults.py
