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
if [ $machine = chloe ]; then
Dirs=(C01 C02 C03 C04 C05 C06 C07 C08 C09 C10)
fi 

for i in ${Dirs[*]}; do	
	echo ' * * * * * * * * * * * * * * * * '$i' * * * * * * * * * * * * * * * * *'
	awk 'NR==26{print;exit}' $i/Out/info.out
	echo '  ========================= run.pipe ==========================='
	tail $n $i/run.pipe
	echo '  ========================= info.out ==========================='
	tail $n $i/Out/info.out
done


