#!/bin/bash
script_path=$(dirname $0)
sim_cwd=$1
sig1=$2
sig2=$3
eps1=$4
eps2=$5
thisName=$6

#echo $script_path
#echo $sim_cwd
#echo $sig1
#echo $sig2
#echo $eps1
#echo $eps2
#echo $thisName

# create working dir and copy all necessary files
mkdir $sim_cwd
sed "s/THISNAMEQQQQ/$thisName/" <$script_path/batch_slurm_prod_1-bromobutane.sh >$sim_cwd/batch_slurm_prod_1-bromobutane.sh 
sed "s/THISNAMEQQQQ/$thisName/" <$script_path/batch_slurm_prod_2-bromobutane.sh >$sim_cwd/batch_slurm_prod_2-bromobutane.sh
cp $script_path/prod_1-bromobutane.mdp $sim_cwd
cp $script_path/prod_2-bromobutane.mdp $sim_cwd
cp $script_path/batch_slurm_eval_PhysProp.sh $sim_cwd

cp $script_path/equilibrated/1-bromobutane_equilibrated.cpt $sim_cwd
cp $script_path/equilibrated/1-bromobutane_equilibrated.gro $sim_cwd
cp $script_path/equilibrated/2-bromobutane_equilibrated.cpt $sim_cwd
cp $script_path/equilibrated/2-bromobutane_equilibrated.gro $sim_cwd
ln -s $script_path/oplsaa_robin.ff $sim_cwd/oplsaa_robin.ff

#adapt parameters
#sed 's/day/night/' <old >new
sed "s/opls_800 opls_800    1/opls_800 opls_800    1    ${sig1:0:7}    ${eps1:0:7}/" <$script_path/equilibrated/1-bromobutane_topol.top >$sim_cwd/topol.tmp
sed "s/opls_730 opls_730    1/opls_730 opls_730    1    ${sig2:0:7}    ${eps2:0:7}/" <$sim_cwd/topol.tmp >$sim_cwd/1-bromobutane_topol.top
rm $sim_cwd/topol.tmp
sed "s/opls_800 opls_800    1/opls_800 opls_800    1    ${sig1:0:7}    ${eps1:0:7}/" <$script_path/equilibrated/2-bromobutane_topol.top >$sim_cwd/topol.tmp
sed "s/opls_730 opls_730    1/opls_730 opls_730    1    ${sig2:0:7}    ${eps2:0:7}/" <$sim_cwd/topol.tmp >$sim_cwd/2-bromobutane_topol.top
rm $sim_cwd/topol.tmp


# execute simulation
cd $sim_cwd
sbatch batch_slurm_prod_1-bromobutane.sh
sleep 2
sbatch batch_slurm_prod_2-bromobutane.sh

