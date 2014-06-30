#!/bin/bash
############################################################################### 

	Dirs=(CDir1 CDir2 CDir3 CDir4 CDir5 CDir6 CDir7 CDir8 CDir9 CDir10 SDir1 SDir2 SDir3 SDir4 SDir5 SDir6 SDir7)

for i in ${Dirs[*]}
do
	
	mv $i/summary.out $i/summary061914.out
	mv $i/InParams.txt $i/InParams061914.txt
	head -1 $i/InParams061914.txt >  $i/InParams.txt
	tail -1 $i/InParams061914.txt >> $i/InParams.txt
#	rm $i/InParams.txt
#	rm $i/summary.out

done

