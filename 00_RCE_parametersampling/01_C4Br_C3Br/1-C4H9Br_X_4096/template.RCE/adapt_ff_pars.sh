script_path=$(dirname $0)
#echo $(pwd) >$script_path/test.txt

sig1=$1
sig2=$2
eps1=$3
eps2=$4


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

sed "s/C800QQQQ  QQQQ/${sig1}e-$sig1exp  ${eps1}e-$eps1exp/" <$script_path/oplsaa_robin.ff/ffnonbonded.itp >$script_path/oplsaa_robin.ff/ffnonbonded.tmp
sed "s/Br730QQQQ  QQQQ/${sig2}e-$sig2exp  ${eps2}e-$eps2exp/" <$script_path/oplsaa_robin.ff/ffnonbonded.tmp >$script_path/oplsaa_robin.ff/ffnonbonded.itp
rm $script_path/oplsaa_robin.ff/ffnonbonded.tmp


