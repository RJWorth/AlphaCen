#!/bin/bash
############################################################################### 
# Run extension $2 on Proxlike/$1 (must first create directory)
# e.g. ./ProxRun.bash Prx01 6-3
# (must be run on the same computer as the original)

pwd=$(pwd)
echo $pwd

prefix='Proxlike/Prx'
#Dirs=(04 05 06 07)			#chloe
Dirs=(02 08 09 10 11)		#shapiro
suffix='/6-3'

for i in ${Dirs[*]}
do
	j=$prefix$i$suffix
	echo $j

	cd $j/Out
		./merc_ext
		./elem
		\mv *.aei AeiOutFiles
	cd $pwd
done


