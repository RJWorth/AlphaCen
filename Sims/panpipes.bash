#!/bin/bash
############################################################################### 
# Write last n lines of each run.pipe file to get most recent changes
# syntax: ./panpipes.sh n d

# choose which directories (with second input d)
if [ $2 = 1 ]; then
	prefix=''
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10 SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 'c' ]; then
	prefix=''
	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10)
fi
if [ $2 = 's' ]; then
	prefix=''
	Dirs=(SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)
fi
if [ $2 = 'p' ]; then
	prefix='Proxlike/Prx'
	Dirs=(01/Disk1-B 02/Disk1-B 03/Disk1-B 04/Disk1-B 05/Disk1-B 01/Disk1-C 02/Disk1-C 03/Disk1-C 04/Disk1-C 05/Disk1-C)
fi
if [ $2 = 'e' ]; then
	prefix='Err/'
	Dirs=(05 07 08 15 19 24 27)
fi

\rm panpipes1.txt

for i in ${Dirs[*]}
do



	echo '******************* '$prefix$i' *******************' >> panpipes1.txt

	if [ $2 = 5 ]; then
		python -c 'import AlphaCenModule as AC; AC.Summary("'$prefix$i'",1e9,wantsum=False,mA=.123,mB=.123)' >> panpipes1.txt
	else
		tail -$1 $prefix$i/run.pipe >> panpipes1.txt
	fi

done

#diff panpipes1.txt panpipes0.txt

