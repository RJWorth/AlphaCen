#!/bin/bash
############################################################################### 
# Predict runtimes for current simulations

if [ $1 = 'C' ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
elif [ $1 = 'S' ]; then
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
else
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7 CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
#	Dirs=(SDir4 CDir3 CDir6)
fi

for i in ${Dirs[*]}
do

	python -c 'import Times;Times.TimeRemaining("'$i'")'

done

#  after above are finished, run finish.sh for summary.out file and R plots
