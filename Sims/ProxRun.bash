#!/bin/bash
############################################################################### 
# Run extension $2 on Proxlike/$1 (must first run ProxSave and edit Dirs below)
# e.g. ./ProxRun.bash
# (must be run on the same computer as the original)

pwd=$(pwd)
echo $pwd

prefix='Err/'	#'Proxlike/Prx'
#Dirs=()		#chloe
Dirs=()		#shapiro
suffix='-copy'

for i in ${Dirs[*]}
do
	j=$prefix$i$suffix
	echo $j

#	gfortran -w -O1 -o $j/Out/merc_$i Files/mercury_TG.for	#j in
#	gfortran -w -O1 -o $j/Out/elem Files/elem.for 			#j in

	cd $j/Out
		./merc_err2
		./elem
		\mv *.aei AeiOutFiles
	cd $pwd
done


