#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=00:20:00
#SBATCH --job-name=1Cemin

gmx grompp -f 01_emin.mdp -c 01_filledbox_1-bromobutane-ep.gro -p topol.top -o 11_emin.tpr
gmx mdrun -deffnm 11_emin
