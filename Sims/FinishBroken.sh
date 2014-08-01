#!/bin/sh
############################################################################### 
### Finish the last piece of a simulation interrupted during the 1e9 stage
### e.g: ./FinishBroken.sh Dir >> Dir/run.pipe &

cd $1/Out
./merc_AC$1
./elem
mv *.aei AeiOutFiles
cd ../..

python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", 1e9, WhichTime="LAST")'

bigstop=$(cat $1/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Proxima-like system.'
		./email.sh $1 'X/X' 'FinishBroken.sh ended - Proxima-like system!'
	else
		./email.sh $1 'X/X' 'FinishBroken.sh ended'
	fi	
