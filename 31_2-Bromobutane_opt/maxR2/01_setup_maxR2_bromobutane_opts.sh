#! /bin/bash

#1
s1="0.220"
s2="0.300"
e1="0.600"
e2="2.500"
wPP="1.2" # sum = 1.2
wQM="0.0 1.2 2.4" # sum = 3.6

#thisCWD="$(pwd)"
#echo $thisCWD
thisCWD="\/home\/rstric2s\/current_sim\/bromobutane.app\/31_2-Bromobutane_opt\/01_C4Br_C3Br\/maxR2"
if [[ -f start_all.sh ]]; then
  rm start_all.sh
fi
touch start_all.sh
chmod +x start_all.sh

if [[ -f AvgRMSEEnergies.csv ]]; then
  rm AvgRMSEEnergies.csv
fi

if [[ -f StatsDensity.csv ]]; then
  rm StatsDensity.csv
fi

if [[ -f StatsDensity_withSimResults.csv ]]; then
  rm StatsDensity_withSimResults.csv
fi

if [[ -f StatsEnergies.csv ]]; then
  rm StatsEnergies.csv
fi

rm slurm-*
rm batch_run_*
rm run_*
rm -rf maxR2-*
rm -rf evalOptParams_DensitySims

for optrun in bestTrainedModels/2-bromo*; do
  optrun="$(cut -d "/" -f2 <<<$optrun | cut -d "_" -f2 | cut -d "." -f1)"
  #echo $optrun
  if [[ -d $optrun ]]; then
    rm -rf $optrun
  fi
  
  cp -r template.maxR2 $optrun
  cd $optrun
  sed "s/cwd: QQQQ/cwd: $thisCWD\/$optrun/" <"bromobutane_hybrid.cfg" >"bromobutane_hybrid.tmp"
  sed "s/logfile: QQQQ/logfile: $optrun.log/" <"bromobutane_hybrid.tmp" >"bromobutane_hybrid.cfg"
  sed "s/weights: QQQQPhysProp/weights: $wPP/" <"bromobutane_hybrid.cfg" >"bromobutane_hybrid.tmp"
  sed "s/weights: QQQQQMMM/weights: $wQM/" <"bromobutane_hybrid.tmp" >"bromobutane_hybrid.cfg"
  rm "bromobutane_hybrid.tmp"
  
  sed "s/SigC800 SigBr730 EpsC800 EpsBr730/$s1 $s2 $e1 $e2/" <parameters.par >parameters.tmp
  mv parameters.tmp parameters.par
  cd ..
  
  #cp bestTrainedModels/1-bromobutane_$optrun.pth $optrun/density/trainedModel_1-bromobutane.pth
  cp bestTrainedModels/2-bromobutane_$optrun.pth $optrun/density/trainedModel_2-bromobutane.pth
  
  n="$(cut -d "-" -f2 <<<$optrun)"
  sed "s/#SBATCH --job-name=2-mR2-QQ/#SBATCH --job-name=2-mR2-$n/" <template.batch_run_maxR2.sh >batch.tmp
  sed "s/.\/run_maxR2-QQ.sh/.\/run_maxR2-$n.sh/" <batch.tmp >batch_run_maxR2-$n.sh
  rm batch.tmp
  echo -e "sbatch batch_run_maxR2-$n.sh" >>start_all.sh
  echo -e "sleep 0.5\n" >>start_all.sh
  
  sed "s/cwd=QQQQ/cwd=$thisCWD\/$optrun/" <template.run_maxR2.sh >run.tmp
  sed "s/logfile=QQQQ/logfile=\"$optrun.log\"/" <run.tmp >run_maxR2-$n.sh
  rm run.tmp
  chmod +x run_maxR2-$n.sh
  #break
done



