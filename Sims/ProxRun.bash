#!/bin/bash
############################################################################### 
# Run extension $2 on Proxlike/$1 (must first run ProxSave and edit Dirs below)
# e.g. ./ProxRun.bash
# (must be run on the same computer as the original)

pwd=$(pwd)
echo $pwd

prefix='Proxlike/Prx'
#Dirs=()			#chloe
Dirs=(12)		#shapiro
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


