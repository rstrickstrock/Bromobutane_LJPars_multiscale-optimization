#! /bin/bash
cwd=QQQQ
logfile=QQQQ

if [[ -d $cwd/PhysProp ]]; then
  rm -rf $cwd/PhysProp
fi
if [[ -d $cwd/QMMM ]]; then
  rm -rf $cwd/QMMM
fi
if [[ -f $cwd/$logfile ]]; then
  rm $cwd/$logfile
fi

miscffoptiw="python /home/rstric2s/current_sim/bromobutane.app/fflow/main.py"

$miscffoptiw $cwd/bromobutane_hybrid.cfg -d
