#!/bin/bash
############################################################################### 

machine=$(hostname -s)
home=$(pwd)
t=0

if [ $machine = Mirka ]; then
Dirs=(01 02)
fi
if [ $machine = shapiro ]; then
Dirs=(01 02 03 04 05 06 07)
fi 

for i in ${Dirs[*]}
do	

	# run sim in 01 directory
	nice -n 10 ./run.bash $i $t > $i/run.pipe &
	echo 'started '$i': '$!

done	# i in Dirs

