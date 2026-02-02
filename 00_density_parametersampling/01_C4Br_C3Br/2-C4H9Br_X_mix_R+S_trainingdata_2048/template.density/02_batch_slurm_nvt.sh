#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=03:00:00
#SBATCH --job-name=1Cnvt

gmx grompp -f 02_nvt.mdp -c 11_emin.gro -p topol.top -o 12_nvt.tpr
gmx mdrun -deffnm 12_nvt
