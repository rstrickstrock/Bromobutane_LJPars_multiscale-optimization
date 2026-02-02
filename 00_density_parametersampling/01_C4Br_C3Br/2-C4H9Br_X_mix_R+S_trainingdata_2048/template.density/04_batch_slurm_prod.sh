#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 50G
#SBATCH --time=72:00:00
#SBATCH --job-name=1Cprod

gmx grompp -f 04_prod.mdp -c 13_npt.gro -t 13_npt.trr -p topol.top -o 14_prod.tpr
gmx mdrun -deffnm 14_prod
