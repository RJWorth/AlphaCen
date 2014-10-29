#!/bin/bash
############################################################################### 
# make a gif of each directory (in series, not parallel)

#Dirs=(Prx12/4-1)
	prefix='Proxlike/Prx'	# e.g. 'Err/' or 'Prx'
#	Dirs=(01 02 03 04 05 06 07 08 09 10 11)
	Dirs=('04/')
	suffix='Disk1-C'

for i in ${Dirs[*]}
do
	### set up directory for this simulation and extension
#	./initialize.bash $prefix$i $suffix

	### run gif-making script for this directory
	./im100.bash $prefix$i$suffix $suffix
	
done

#./GifEmail.sh $i 'gif'

