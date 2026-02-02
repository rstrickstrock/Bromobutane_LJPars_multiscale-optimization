#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalRCE


cd output/

if [ -f Energy.abs.kcal.txt ]; then
  rm Energy.abs.kcal.txt
fi
touch Energy.abs.kcal.txt

pot=$(echo "Potential" | gmx energy -f 21_emin.edr -o 21_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "2-bromobutane_conformer_01 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 22_emin.edr -o 22_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "2-bromobutane_conformer_02 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 23_emin.edr -o 23_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "2-bromobutane_conformer_03 $pot" >>Energy.abs.kcal.txt
sleep 1

cd ..
python calcRCE.py
