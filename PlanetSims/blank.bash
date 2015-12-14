#!/bin/bash
###################################################################################

machine=$(hostname -s)

if [ $machine = Mirka ]; then
Dirs=(C41)
elif [ $machine = chloe ];then
Dirs=()
elif [ $machine = shapiro ]; then
Dirs=()
fi

for i in ${Dirs[*]}
do	

	echo $i

	mkdir $i
	mkdir $i/In
	mkdir $i/Out
	mkdir $i/Aei

	./pipe.bash $i

done






