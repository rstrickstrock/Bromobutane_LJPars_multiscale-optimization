#! /bin/bash

rndInts=("678 147" "561 237" "588 951" "490 395" "877 297" "721 711" "985 171" "75 16" "669 530" "999 794" "936 111" "816 968" "48 986" "829 996" "272 759" "390 930" "633 928" "854 554" "562 78" "222 294" "725 582" "731 249" "791 35" "180 510" "593 634")

if [[ -f 01_run-all.sh ]]; then
  rm 01_run-all.sh
fi

touch 01_run-all.sh

if [[ -f 01_batch_NNTraining-* ]]; then
  rm 01_batch_NNTraining-*
fi

if [[ -f 01_NNTraining-* ]]; then
  rm 01_NNTraining-*
fi

if [[ -f slurm-* ]]; then
  rm slurm-*
fi

if [[ -d trainedModels ]]; then
  rm -rf trainedModels
fi

i=0
for rndInt in "${rndInts[@]}"; do
  #echo $rndInt
  SPLITEDrndInt=(${rndInt//' '/ })
  A=${SPLITEDrndInt[0]}
  B=${SPLITEDrndInt[1]}
  #echo $FIRSTRANDOMINTEGER
  #echo $SECONDRANDOMINTEGER
  if [ ${#i} -eq 1 ]; then
    k="0$i"
  else
    k=$i
  fi 
  if [[ -f 01_NNTraining-$k.py ]]; then
    rm 01_NNTraining-$k.py
  fi
  if [[ -f 01_batch_NNTraining-$k.sh ]]; then
    rm 01_batch_NNTraining-$k.sh
  fi
  sed "s/rndInts = \[miauA, miauB\]/rndInts = \[$A, $B\]/g" <00_NNTraining-X.py >01_NNTraining-$k.py
  sed "s/job-name=NN-X/job-name=NN-$k/g" <00_batch_NNTraining-X.sh >tmpBatchFile
  sed "s/01_NNTraining-X.py/01_NNTraining-$k.py/g" <tmpBatchFile >01_batch_NNTraining-$k.sh
  rm tmpBatchFile
  echo "sbatch 01_batch_NNTraining-$k.sh" >>01_run-all.sh
  echo "sleep 2" >>01_run-all.sh
  let "i++"
done
chmod +x 01_run-all.sh
