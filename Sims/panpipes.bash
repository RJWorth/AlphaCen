#!/bin/bash
############################################################################### 
# Write last few lines of each run.pipe file to get most recent changes
# syntax: ./panpipes.sh n d

# choose which directories (with second input d)
if   [ $2 = 1 ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 SDir1 SDir2 SDir3 CDir10 NDir1)
fi
if [ $2 = 2 ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7)
fi
if [ $2 = 3 ]; then
	Dirs=(SDir1 SDir2 SDir3)
fi
if [ $2 = 3 ]; then
	Dirs=(NDir1 CDir10)
fi


# write last n lines, n = first input
for i in ${Dirs[*]}
do
	echo $i
	tail -$1 $i/run.pipe
done

