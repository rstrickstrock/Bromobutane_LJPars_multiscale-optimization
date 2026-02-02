#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --time=02:00:00
#SBATCH --job-name=1THISNAMEQQQQ

gmx grompp -f prod_1-bromobutane.mdp -c 1-bromobutane_equilibrated.gro -t 1-bromobutane_equilibrated.cpt -p 1-bromobutane_topol.top -o 1-bromobutane_prod.tpr

gmx mdrun -v -deffnm 1-bromobutane_prod

