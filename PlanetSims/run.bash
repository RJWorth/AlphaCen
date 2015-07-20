#!/bin/bash
############################################################################### 

home=$(pwd)
Dirs=(01)


for i in ${Dirs[*]}
do
	# write param.in file
	python -c "import Merc as M; M.WriteParamInFile(loc='"$i"/In/', tf=365.25e0, mStar=1.105, rEj=100, rStar=0.005)"

	# write big.in file with stars+planetesimals
	# (currently only has planetesimals around AlCenA)
	python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(), loc='"$i"/')"

	# run sim in 01 directory
	cd $i
	./run.sh
	cd ..
done	# i in Dirs

