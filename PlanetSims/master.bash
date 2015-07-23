#!/bin/bash
############################################################################### 

machine=$(hostname -s)
home=$(pwd)
ti=0
tf=0

if [ $machine = Mirka ]; then
Dirs=(L01)
fi
if [ $machine = shapiro ]; then
Dirs=(S01 S02 S03 S04 S05 S06 S07)
fi 

for i in ${Dirs[*]}
do	

	# run sim in 01 directory
	nice -n 10 ./run.bash $i $ti $tf > $i/run.pipe &
	echo 'started '$i': '$!

done	# i in Dirs

