#!/bin/bash
############################################################################### 
wd=$(pwd)
echo $wd

echo 'do100, dir = '$1
suffix=${1:(-3)}

### Generate images in R
echo $1', reading AEI'
R CMD BATCH -$1 ReadAei.R
	mv gif/$1/CoordvsT.pdf gif/$1/CoordvsT$suffix.pdf 
	mv gif/$1/XY.pdf gif/$1/XY$suffix.pdf 

### Read in number of timesteps
lengthT=$(cat lengthT.txt)
step1=$(head -1 steps.txt)
steps=$(tail -1 steps.txt)

### Write the actual script for each version
echo $1', writing doim.sh'
for i in 2 3; do
### first line
	echo $(pwd)' '$i
	name=$(printf "%04d" $step1)
	echo 'convert -delay 10  -size 560x560 gifimgs/'$i'xy'$name'.jpg \' > img/$1/doim.sh
### all subsequent lines
	#for ((x=2;x<=$lengthT;x+=1)); 
	for x in $steps; do 
		name=$(printf "%04d" $x)
		echo '                     -page +0+0  gifimgs/'$i'xy'$name'.jpg \' >> img/$1/doim.sh
		done
	echo '         -loop 0  '$wd'/gif/'$1'/'$i'xy'$suffix'.gif' >> img/$1/doim.sh
	chmod a+x img/$1/doim.sh

### Run that script
	echo $1', running convert'
	echo $(pwd)' before doim.sh '$i
	cd img/$1; ./doim.sh; cd $wd
	echo $(pwd)'  after doim.sh '$i
#		if [[ $1 == */* ]]; then
#		cd ..;
#		fi
done

#echo $1', moving gifs'
#\cp -p gif/xy_$1.gif ~/AlphaCen/Analysis/gif
#eog gif/ae_$1.gif



