#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:20:00
#SBATCH --job-name=RCE

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_01.gro -p topol_1-bromobutane-X.top -o output/11_nvt.tpr
gmx mdrun -deffnm output/11_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_02.gro -p topol_1-bromobutane-X.top -o output/12_nvt.tpr
gmx mdrun -deffnm output/12_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_03.gro -p topol_1-bromobutane-X.top -o output/13_nvt.tpr
gmx mdrun -deffnm output/13_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_04.gro -p topol_1-bromobutane-X.top -o output/14_nvt.tpr
gmx mdrun -deffnm output/14_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_05.gro -p topol_1-bromobutane-X.top -o output/15_nvt.tpr
gmx mdrun -deffnm output/15_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_06.gro -p topol_1-bromobutane-X.top -o output/16_nvt.tpr
gmx mdrun -deffnm output/16_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_07.gro -p topol_1-bromobutane-X.top -o output/17_nvt.tpr
gmx mdrun -deffnm output/17_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_08.gro -p topol_1-bromobutane-X.top -o output/18_nvt.tpr
gmx mdrun -deffnm output/18_nvt
sleep 2

gmx grompp -f 00_nvt.mdp -c 1-bromobutane-X_conformer_09.gro -p topol_1-bromobutane-X.top -o output/19_nvt.tpr
gmx mdrun -deffnm output/19_nvt

