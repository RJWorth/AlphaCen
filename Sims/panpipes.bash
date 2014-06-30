#!/bin/bash
############################################################################### 
# Write last few lines of each run.pipe file to get most recent changes
# syntax: ./panpipes.sh n d

# choose which directories (with second input d)
if [ $2 = 1 ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10 SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 2 ]; then
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
fi
if [ $2 = 3 ]; then
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 4 ]; then
	Dirs=(PDir1 PDir2 PDir3)
fi


# write last n lines, n = first input

\rm panpipes1.txt

for i in ${Dirs[*]}
do
	echo $i >> panpipes1.txt
	tail -$1 $i/run.pipe >> panpipes1.txt
done

#diff panpipes1.txt panpipes0.txt

