#!/bin/bash
##############################################################################

machine=rjw274@chloe.astro.psu.edu
dir=C09

scp -rp $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Aei' $dir
scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/In/big.in' $dir/In
scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Out/info.out' $dir/Out
scp -p  $machine:'/Volumes/Macintosh\ HD\ 2/rjw/ACPltSims/'$dir'/Out/element.out' $dir/Out


