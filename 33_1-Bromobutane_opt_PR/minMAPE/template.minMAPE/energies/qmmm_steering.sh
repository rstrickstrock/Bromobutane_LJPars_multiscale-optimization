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

sig1exp=0
while [ $(echo "$sig1 < 1.0" | bc ) -eq 1 ]; do
  sig1=$(echo $sig1*10 | bc )
  sig1exp=$((sig1exp+1))
done
sig1="${sig1}0000000"
sig1=${sig1:0:7}
if [[ ${#sig1exp} -lt 10 ]]; then
  sig1exp="0$sig1exp"
fi

sig2exp=0
while [ $(echo "$sig2 < 1.0" | bc ) -eq 1 ]; do
  sig2=$(echo $sig2*10 | bc )
  sig2exp=$((sig2exp+1))
done
sig2="${sig2}0000000"
sig2=${sig2:0:7}
if [[ ${#sig2exp} -lt 10 ]]; then
  sig2exp="0$sig2exp"
fi

eps1exp=0
while [ $(echo "$eps1 < 1.0" | bc ) -eq 1 ]; do
  eps1=$(echo $eps1*10 | bc )
  eps1exp=$((eps1exp+1))
done
eps1="${eps1}0000000"
eps1=${eps1:0:7}
if [[ ${#eps1exp} -lt 10 ]]; then
  eps1exp="0$eps1exp"
fi

eps2exp=0
while [ $(echo "$eps2 < 1.0" | bc ) -eq 1 ]; do
  eps2=$(echo $eps2*10 | bc )
  eps2exp=$((eps2exp+1))
done
eps2="${eps2}0000000"
eps2=${eps2:0:7}
if [[ ${#eps2exp} -lt 10 ]]; then
  eps2exp="0$eps2exp"
fi

sim_cwd=$sim_cwd/$sim_type.$iteration.$direction
# create working dir and copy all necessary files
mkdir $sim_cwd
mkdir $sim_cwd/output/

cp $script_path/batch_slurm_RCE.sh $sim_cwd/
cp $script_path/00_emin.mdp $sim_cwd/
cp $script_path/1-bromobutane-X_conformer_0* $sim_cwd/
#cp $script_path/2-bromobutane-X_conformer_0* $sim_cwd/
cp $script_path/topol_1-bromobutane-X.top $sim_cwd/
#cp $script_path/topol_2-bromobutane-X.top $sim_cwd/
cp -r $script_path/oplsaa_robin.ff/ $sim_cwd/

#adapt parameters
#sed "s/opls_136 opls_136    1/opls_136 opls_136    1    ${sig1:0:7}    ${eps1:0:7}/" <$script_path/topol.top >$sim_cwd/topol.tmp
#sed "s/opls_730 opls_730    1/opls_730 opls_730    1    ${sig2:0:7}    ${eps2:0:7}/" <$sim_cwd/topol.tmp >$sim_cwd/topol.top
#rm $sim_cwd/topol.tmp

sed "s/C800QQQQ  QQQQ/${sig1}e-$sig1exp  ${eps1}e-$eps1exp/" <$sim_cwd/oplsaa_robin.ff/ffnonbonded.itp >$sim_cwd/oplsaa_robin.ff/ffnonbonded.tmp
sed "s/Br730QQQQ  QQQQ/${sig2}e-$sig2exp  ${eps2}e-$eps2exp/" <$sim_cwd/oplsaa_robin.ff/ffnonbonded.tmp >$sim_cwd/oplsaa_robin.ff/ffnonbonded.itp
rm $sim_cwd/oplsaa_robin.ff/ffnonbonded.tmp

#copy files for simulation evaluation
cp $script_path/batch_slurm_eval_RCE.sh $sim_cwd
cp $script_path/calcRCE.py $sim_cwd

# execute simulation
cd $sim_cwd
sbatch batch_slurm_RCE.sh
