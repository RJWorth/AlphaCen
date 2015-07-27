#!/bin/bash
############################################################################### 

machine=$(hostname -s)
home=$(pwd)

### Sim parameters
ti=0
tf=6
rtr=3.08
rice=2.7
alpha=1.5
sigma=3.

if [ $machine = Mirka ]; then
Dirs=(L01)
fi
if [ $machine = shapiro ]; then
Dirs=(S01 S02 S03 S04 S05 S06 S07)
fi
if [ $machine = chloe ];then
Dirs=(C01)
fi 

for i in ${Dirs[*]}
do	

	# run sim in 01 directory
	nice -n 10 ./run.bash $i $ti $tf $rtr $sigma > $i/run.pipe 2>&1 &
	echo 'started '$i': '$!

done	# i in Dirs

