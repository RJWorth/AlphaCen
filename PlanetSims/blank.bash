#!/bin/bash
###################################################################################

machine=$(hostname -s)

if [ $machine = Mirka ]; then
Dirs=(Bin25 Bin26 Bin27 Bin28 Bin29 Bin30 Bin31 Bin32 Bin33 Bin34 Bin35 Bin36)
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






