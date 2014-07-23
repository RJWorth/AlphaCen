#!/bin/bash
############################################################################### 

	prefix='Proxlike/Prx'
	Dirs=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22)
#	Dirs=(01)
#	suffix='/6-3'


for i in ${Dirs[*]}
do
	j=$prefix$i$suffix
	echo $j

### Generate summary/plots
	python -c 'import AlphaCenModule;AlphaCenModule.Summary("'$j'",1e9,wantplot=True,wantsum=False)'

done

