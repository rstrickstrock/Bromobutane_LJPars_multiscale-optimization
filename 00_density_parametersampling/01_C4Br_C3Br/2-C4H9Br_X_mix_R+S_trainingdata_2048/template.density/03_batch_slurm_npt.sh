#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 10G
#SBATCH --time=10:00:00
#SBATCH --job-name=1Cnpt

gmx grompp -f 03_npt.mdp -c 12_nvt.gro -t 12_nvt.trr -p topol.top -o 13_npt.tpr
gmx mdrun -deffnm 13_npt
