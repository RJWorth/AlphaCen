#!/bin/bash
###############################################################################
# scp all the specified directories into the analogous location here

pre='Proxlike/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
Dir=(11 12 13 14)
suf=('/DiskB-2' '/DiskB-3')

for i in ${Dir[*]}; do
	echo $pre$i
	scp -p rjw274@hammer.rcc.psu.edu:'~/work/AlphaCen/Sims/'$pre$i'/DiskSurv.pdf' $pre$i
	scp -p rjw274@hammer.rcc.psu.edu:'~/work/AlphaCen/Sims/'$pre$i'/SurvGrid.txt' $pre$i

	for j in ${suf[*]}; do
		d=$pre$i$j
#		mkdir -p $d
		echo $d
		scp -p rjw274@hammer.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d'/DiskSurv.pdf' $d

#		echo R CMD BATCH -$d '../Analysis/ReadDisk.R'
#		R CMD BATCH -$d '../Analysis/ReadDisk.R'

	done
done
