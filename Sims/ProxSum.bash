#!/bin/bash
############################################################################### 
# run Summary on all Prx, and put plots in one directory
# e.g: ./ProxSum.bash

	prefix='Proxlike/Prx'
	Dirs=(01 02 03 04 05 06 07 08 09 10 11)
#	Dirs=(01 02)
	suffix='/Original'

\rm Proxlike/Plots/sumfeed.txt
for i in ${Dirs[*]}
do
	j=$prefix$i$suffix
	echo $j

### Generate summary/plots
	echo '*********************'$j'*********************' >> Proxlike/Plots/sumfeed.txt
	python -c 'import AlphaCenModule;AlphaCenModule.Summary("'$j'",1e9,WhichTime="'$i'",wantplot=True,wantsum=True)' >> Proxlike/Plots/sumfeed.txt
	mv $j/EvsT_AB.png  Proxlike/Plots/EvsT_AB_$i.png
	mv $j/EvsT_ABC.png Proxlike/Plots/EvsT_ABC_$i.png
	mv $j/RvsT_AB.png  Proxlike/Plots/RvsT_AB_$i.png
	mv $j/RvsT_ABC.png Proxlike/Plots/RvsT_ABC_$i.png
	mv $j/EvsR_AB.png  Proxlike/Plots/EvsR_AB_$i.png
	mv $j/EvsR_ABC.png Proxlike/Plots/EvsR_ABC_$i.png
done

