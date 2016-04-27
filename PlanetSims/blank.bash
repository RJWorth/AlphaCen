#!/bin/bash
###################################################################################

machine=$(hostname -s)

if [ $machine = Mirka ]; then
Dirs=(Bin14 Bin15 Bin16 Bin17 Bin18 Bin19 Bin20 Bin21 Bin22 Bin23 Bin24)
elif [ $machine = chloe ];then
Dirs=()
elif [ $machine = shapiro ]; then
Dirs=()
fi

for i in ${Dirs[*]}
do	

	echo $i

	./pipe.bash $i

done






