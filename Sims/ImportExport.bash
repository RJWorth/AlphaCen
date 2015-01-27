#!/bin/bash
###############################################################################
# compare the orbital parameters in the original and disk triple systems

#module load python/2.7.3
#module load gcc/4.7.1
#module load R

host=rjw274@hammer.rcc.psu.edu

h=$(pwd)
echo $h

### Older true-mass batches of sims
pre='Proxlike/071714/Prx'
Dir=(1 2 3 4 5 6 7 8 9 10 11 12)
suf='-True'
tag='0717'


#pre='Proxlike/072314/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)
#suf='-True'
#tag='0723'

#pre='Proxlike/081414/Prx'
#Dir=(01 02 03 04 05 06 07 08)
#suf='-True'
#tag='0814'

### More recent equal-mass batch
#pre='Proxlike/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
#Dir=(01 03 05 06 07 09 10 11 12 13 14 15 16 17 19 20 21 22 23 24 26 33)
#suf='-Equal'
#tag=''

#for i in ${Dir[*]}; do

#	d=$pre$i
#	echo '--------------------------------'$d'--------------------------------'
#	scp -p $host:'~/work/AlphaCen/Sims/'$d'/DiskSurv.eps' Proxlike/Plots/DiskSurv$suf$tag'-'$i.eps

#done

#gs -dBATCH -dNOPAUSE -dPDFFitPage -dSAFER -q -dDEVICEWIDTHPOINTS=360 -dDEVICEHEIGHTPOINTS=270 -sDEVICE=pdfwrite -sOutputFile=Proxlike/Plots/All-DiskSurv${suf}.pdf Proxlike/Plots/DiskSurv${suf}*.eps

### Update all DiskSummary files
#if [ $suf = '-Equal' ]; then
	scp -p $host:'~/work/AlphaCen/Sims/Proxlike/DiskSummary.txt' Proxlike/Plots/DiskSummary-Equal.txt
	scp -p $host:'~/work/AlphaCen/Sims/Proxlike/PrxDisksSurvival.pdf' Proxlike/Plots/PrxDisksSurvival.eps
#elif [ $tag = '0717' ]; then
	scp -p $host:'~/work/AlphaCen/Sims/Proxlike/071714/DiskSummary.txt' Proxlike/Plots/DiskSummary-0717.txt
#elif [ $tag = '0723' ]; then
	scp -p $host:'~/work/AlphaCen/Sims/Proxlike/072314/DiskSummary.txt' Proxlike/Plots/DiskSummary-0723.txt
#elif [ $tag = '0814' ]; then
	scp -p $host:'~/work/AlphaCen/Sims/Proxlike/081414/DiskSummary.txt' Proxlike/Plots/DiskSummary-0814.txt
#fi

