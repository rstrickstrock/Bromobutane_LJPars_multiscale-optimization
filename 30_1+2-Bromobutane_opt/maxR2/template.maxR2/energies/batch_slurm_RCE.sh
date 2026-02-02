#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 1G
#SBATCH --time=00:20:00
#SBATCH --job-name=RCE

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_01.gro -p topol_1-bromobutane-X.top -o output/11_emin.tpr
gmx mdrun -deffnm output/11_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_02.gro -p topol_1-bromobutane-X.top -o output/12_emin.tpr
gmx mdrun -deffnm output/12_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_03.gro -p topol_1-bromobutane-X.top -o output/13_emin.tpr
gmx mdrun -deffnm output/13_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_04.gro -p topol_1-bromobutane-X.top -o output/14_emin.tpr
gmx mdrun -deffnm output/14_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_05.gro -p topol_1-bromobutane-X.top -o output/15_emin.tpr
gmx mdrun -deffnm output/15_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_06.gro -p topol_1-bromobutane-X.top -o output/16_emin.tpr
gmx mdrun -deffnm output/16_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_07.gro -p topol_1-bromobutane-X.top -o output/17_emin.tpr
gmx mdrun -deffnm output/17_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_08.gro -p topol_1-bromobutane-X.top -o output/18_emin.tpr
gmx mdrun -deffnm output/18_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 1-bromobutane-X_conformer_09.gro -p topol_1-bromobutane-X.top -o output/19_emin.tpr
gmx mdrun -deffnm output/19_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_01.gro -p topol_2-bromobutane-X.top -o output/21_emin.tpr
gmx mdrun -deffnm output/21_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_02.gro -p topol_2-bromobutane-X.top -o output/22_emin.tpr
gmx mdrun -deffnm output/22_emin
sleep 2

gmx grompp -f 00_emin.mdp -c 2-bromobutane-X_conformer_03.gro -p topol_2-bromobutane-X.top -o output/23_emin.tpr
gmx mdrun -deffnm output/23_emin
