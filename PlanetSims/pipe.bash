#!/bin/bash
##############################################################################

machine=rjw274@chloe.astro.psu.edu
dir=Data/$1

### make directories for this sim
#mkdir $dir
#mkdir $dir/In
#mkdir $dir/Out

### download individual files for each sim, if not already present
files=(In/big.in In/small.in In/param.in Out/info.out Out/element.out Out/param.dmp)
for f in ${files[*]}; do
	if [ ! -f $dir/$f ];then
#	scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/'$f $dir/$f
	echo $dir'/'$f 'missing'
	fi
done

### download the .aei directory, if not already present
if [ ! -f $dir/Aei/AlCenB.aei ];then
#	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei' $dir
	echo $dir/aei missing
fi

