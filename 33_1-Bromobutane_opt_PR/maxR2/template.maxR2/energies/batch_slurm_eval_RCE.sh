#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:10:00
#SBATCH --job-name=evalRCE

cd output

if [ -f Energy.abs.kcal.txt ]; then
  rm Energy.abs.kcal.txt
fi
touch Energy.abs.kcal.txt

pot=$(echo "Potential" | gmx energy -f 11_emin.edr -o 11_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_01 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 12_emin.edr -o 12_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_02 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 13_emin.edr -o 13_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_03 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 14_emin.edr -o 14_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_04 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 15_emin.edr -o 15_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_05 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 16_emin.edr -o 16_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_06 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 17_emin.edr -o 17_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_07 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 18_emin.edr -o 18_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_08 $pot" >>Energy.abs.kcal.txt
sleep 1

pot=$(echo "Potential" | gmx energy -f 19_emin.edr -o 19_potential | grep "Potential")
sleep 1
pot=$(echo $pot | cut -f2 -d" ")
echo -e "1-bromobutane_conformer_09 $pot" >>Energy.abs.kcal.txt
sleep 1

cd ..
python calcRCE.py
