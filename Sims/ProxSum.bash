#!/bin/bash
############################################################################### 
# run Summary on all Prx, and put plots in one directory
# e.g: ./ProxSum.bash

	prefix='Proxlike/080614/Prx'
	Dirs=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16)
#	Dirs=(02)
	suffix='/Original'

#\rm Proxlike/Plots/sumfeed.txt
for i in ${Dirs[*]}
do
	j=$prefix$i$suffix
	echo $j

### Generate summary/plots
#	echo '*********************'$j'*********************' >> Proxlike/Plots/sumfeed.txt
#	echo ' - running Summary'
#	python -c 'import AlphaCenModule;AlphaCenModule.Summary("'$j'",1e9,WhichTime="'$i'",wantplot=True,wantsum=True)' >> Proxlike/Plots/sumfeed.txt

### Generate aei file
#	echo ' - running WriteAEI'
#	python -c 'import AlphaCenModule;AlphaCenModule.WriteAEI("'$j'",1e9)'

### Make R pdf plots from AEI files
	echo ' - making R plots'
	R CMD BATCH -$j TimeData.R

### Convert R plots from pdf to jpg, put in one directory
	echo ' - converting plots'
	convert $j/Out/AeiOutFiles/EpsAEIvsT.pdf  Proxlike/Plots/EpsAEIvsT_$i.jpg
	convert $j/Out/AeiOutFiles/EvsA.pdf       Proxlike/Plots/EvsA_$i.jpg
done

