#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalRCE

echo "Potential" | gmx energy -f output/11_emin.edr -o output/11_energy >>output/11_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/12_emin.edr -o output/12_energy >>output/12_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/13_emin.edr -o output/13_energy >>output/13_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/14_emin.edr -o output/14_energy >>output/14_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/15_emin.edr -o output/15_energy >>output/15_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/16_emin.edr -o output/16_energy >>output/16_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/17_emin.edr -o output/17_energy >>output/17_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/18_emin.edr -o output/18_energy >>output/18_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/19_emin.edr -o output/19_energy >>output/19_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/21_emin.edr -o output/21_energy >>output/21_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/22_emin.edr -o output/22_energy >>output/22_potentialEnergy.txt
sleep 2

echo "Potential" | gmx energy -f output/23_emin.edr -o output/23_energy >>output/23_potentialEnergy.txt
sleep 2

python gatherResults.py
