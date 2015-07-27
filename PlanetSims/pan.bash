#!/bin/bash
##############################################################################

machine=$(hostname -s)
n=$1

if [ $machine = Mirka ]; then
Dirs=(L01)
fi
if [ $machine = shapiro ]; then
Dirs=(S01 S02 S03 S04 S05 S06 S07)
fi 

for i in ${Dirs[*]}; do	
	echo '=========================== '$i' ============================'
	echo '------------------------- run.pipe --------------------------'
	tail $n $i/run.pipe
	echo '------------------------- info.out --------------------------'
	tail $n $i/Out/info.out
done


