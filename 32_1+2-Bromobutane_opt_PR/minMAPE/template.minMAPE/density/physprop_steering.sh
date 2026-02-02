#!/bin/bash
script_path=$(dirname $0)
sim_cwd=$1
sim_type=$2
iteration=$3
direction=$4
sig1=$5
sig2=$6
eps1=$7
eps2=$8

sim_cwd=$sim_cwd/$sim_type.$iteration.$direction
#echo $sim_cwd
# create working dir and copy all necessary files
mkdir $sim_cwd
cp $script_path/batch_slurm_evalSurrogateModel.sh $sim_cwd
cp $script_path/loadEvalSurrogModel.py $sim_cwd
cp $script_path/trainedModel_1-bromobutane.sav $sim_cwd
cp $script_path/trainedModel_2-bromobutane.sav $sim_cwd
cp $script_path/modelConfig.csv $sim_cwd

#adapt parameters
#sed 's/day/night/' <old >new
sed "s/SC800,SBr730,EC800,EBr730/$sig1,$sig2,$eps1,$eps2/" <$script_path/this_parameters.csv >$sim_cwd/this_parameters.csv

#copy files for simulation evaluation
cp $script_path/batch_slurm_eval_PhysProp.sh $sim_cwd

# execute simulation
cd $sim_cwd
sbatch batch_slurm_evalSurrogateModel.sh

