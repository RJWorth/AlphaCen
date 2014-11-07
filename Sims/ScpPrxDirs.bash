#!/bin/bash
###############################################################################
# scp all the specified directories into the analogous location here

pre='Proxlike/Prx'
Dir=(01) # 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
suf=('/DiskA-2' '/DiskA-3' '/DiskB-2' '/DiskB-3')

for i in ${Dir[*]}; do
	for j in ${suf[*]}; do
		d=$pre$i$j
		mkdir -p $d
		echo    rjw274@lionxj.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d/*.pdf $d
		scp -rp rjw274@lionxj.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d $d 		
#		scp -rp rjw274@shapiro.astro.psu.edu:'~/AlphaCen/Sims/'$d $pre$i
	done
done