#!/bin/bash
###################################################################################

machine=$(hostname -s)

if [ $machine = Mirka ]; then
Dirs=(C01 C05 C06 C09)
Dirs=(C02 C03 C04 C07 C08 C10 C11 C12)
elif [ $machine = chloe ];then
#Dirs=(C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12)
Dirs=(C03 C04 C07 C08 C10 C11 C12)
elif [ $machine = shapiro ]; then
Dirs=(S01 S02 S03 S04 S05 S06 S07)
fi

for i in ${Dirs[*]}
do	

	echo $i

#	mkdir $i
#	mkdir $i/In
#	mkdir $i/Out
#	mkdir $i/Aei

	./pipe.bash $i

done






