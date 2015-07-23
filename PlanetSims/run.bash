#!/bin/bash
############################################################################### 

home=$(pwd)
Dirs=(01 02 03 04 05 06 07)
t=2

for i in ${Dirs[*]}
do	
	loc=$i/In/
	# write param.in file
	python -c "import Merc as M; M.WriteParamInFile(loc='"$loc"', tf=365.25e'"$t"', mStar=1.105, rEj=100, rStar=0.005)"

	# write big.in file with stars+planetesimals
	# (currently only has planetesimals around AlCenA)
	python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(), loc='"$loc"')"

	# run sim in 01 directory
	cd $i
	./run.sh > run.pipe &
	echo 'started '$i': '$!
	cd ..
done	# i in Dirs

