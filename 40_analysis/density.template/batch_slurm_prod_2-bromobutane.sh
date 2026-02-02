#!/bin/bash
#SBATCH --partition=hpc,hpc1,hpc3
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --time=02:00:00
#SBATCH --job-name=2THISNAMEQQQQ

gmx grompp -f prod_2-bromobutane.mdp -c 2-bromobutane_equilibrated.gro -t 2-bromobutane_equilibrated.cpt -p 2-bromobutane_topol.top -o 2-bromobutane_prod.tpr

gmx mdrun -v -deffnm 2-bromobutane_prod

