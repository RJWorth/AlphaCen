#!/bin/bash
##############################################################################

machine=rjw274@chloe.astro.psu.edu
dir=$1

scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/In/big.in' $dir/In
scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Out/info.out' $dir/Out
scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Out/element.out' $dir/Out

### download all, or a subset, of the .aei files

if [ $dir = 'C01' ];then
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/AlCenA.aei' $dir/Aei
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/P0018.aei' $dir/Aei
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/P0053.aei' $dir/Aei
elif [ $dir = 'C06' ];then
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/AlCenA.aei' $dir/Aei
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/P0003.aei' $dir/Aei
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei/P0911.aei' $dir/Aei
else
	scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei' $dir
fi
