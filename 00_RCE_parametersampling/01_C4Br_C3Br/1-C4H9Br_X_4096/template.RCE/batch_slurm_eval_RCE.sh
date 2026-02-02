#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalRCE

echo "Total-Energy" | gmx energy -f output/11_nvt.edr -o output/11_energy >>output/11_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/12_nvt.edr -o output/12_energy >>output/12_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/13_nvt.edr -o output/13_energy >>output/13_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/14_nvt.edr -o output/14_energy >>output/14_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/15_nvt.edr -o output/15_energy >>output/15_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/16_nvt.edr -o output/16_energy >>output/16_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/17_nvt.edr -o output/17_energy >>output/17_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/18_nvt.edr -o output/18_energy >>output/18_totalEnergy.txt
sleep 2

echo "Total-Energy" | gmx energy -f output/19_nvt.edr -o output/19_energy >>output/19_totalEnergy.txt
sleep 2

python gatherResults.py
