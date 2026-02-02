#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:20:00
#SBATCH --job-name=RCE


gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_01.gro -p topol_2-bromobutane-X.top -o output/21_emin.tpr
gmx mdrun -deffnm output/21_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_02.gro -p topol_2-bromobutane-X.top -o output/22_emin.tpr
gmx mdrun -deffnm output/22_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_03.gro -p topol_2-bromobutane-X.top -o output/23_emin.tpr
gmx mdrun -deffnm output/23_emin
