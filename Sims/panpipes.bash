#!/bin/bash
############################################################################### 
# Write last few lines of each run.pipe file to get most recent changes
# syntax: ./panpipes.sh n d

# choose which directories (with second input d)
if [ $2 = 1 ]; then
	prefix=''
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10 SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 2 ]; then
	prefix=''
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
fi
if [ $2 = 3 ]; then
	prefix=''
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 4 ]; then
	prefix='Proxlike/Prx'
	Dirs=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22)
fi


# write last n lines, n = first input

\rm panpipes1.txt

for i in ${Dirs[*]}
do
	echo $prefix$i >> panpipes1.txt
	tail -$1 $prefix$i/run.pipe >> panpipes1.txt
done

#diff panpipes1.txt panpipes0.txt

