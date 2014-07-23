#!/bin/bash
############################################################################### 
# make a gif of each directory (in series, not parallel)

#Dirs=(Prx12/4-1)
	prefix='Proxlike/Prx'
#	Dirs=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22)
	Dirs=(01)
	suffix='6-3'

for i in ${Dirs[*]}
do
	# run gif-making script for this directory
	./im100.bash $prefix$i/$suffix
	
done

#./GifEmail.sh $i 'gif'

