#!/bin/bash
############################################################################### 
# make a gif of each directory (in series, not parallel)

#Dirs=(7)
Dirs=(Prx12/4-1)

for i in ${Dirs[*]}
do
	# if first time:
#	mkdir img/$Dirs
#	mkdir img/$Dirs/gifimgs	
#	mkdir gif/$Dirs

	# run gif-making script for this directory
	./im100.bash $i
	
done

#./GifEmail.sh $i 'gif'

