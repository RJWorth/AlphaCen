#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir

# to clean previous run:
# move plots, SumAll, AllParams, Rout files to Saved/whatever
# \cp -p BlankDir/summary.out *Dir*		or
# \rm *Dir*/summary.out
# \rm *Dir*/InParams.txt

machine=$(hostname -s)
if [ $machine = chloe ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
	Dirs=()
fi
if [ $machine = shapiro ]; then
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
	Dirs=()
fi

for i in ${Dirs[*]}
do

#	nice -n 10 ./FinishBroken.sh $i 'equal' >> $i/run.pipe &
	
	nice -n 10 ./run.sh $i > $i/run.pipe &
	echo 'master: '$i'  '$!

done

#  after above are finished, run finish.sh for SumAll.out file and R plots
