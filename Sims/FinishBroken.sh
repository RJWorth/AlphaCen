#!/bin/sh
############################################################################### 
### Finish the last piece of a simulation interrupted during the 1e9 stage
### e.g: ./FinishBroken.sh Dir flag >> Dir/run.pipe &

cd $1/Out
./merc_AC$1
./elem
mv *.aei AeiOutFiles
cd ../..

### Run summary (pick appropriate masses)
if [ $2 = 'unequal' ]; then
	python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", 1e9, WhichTime="LAST", mA=1.105, mB=.934, mC=.123)'
elif [ $2 = 'equal' ]; then
	python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", 1e9, WhichTime="LAST", mA=.123, mB=.123, mC=.123)'
else
	echo 'invalid mass flag'

bigstop=$(cat $1/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Proxima-like system.'
		./email.sh $1 'X/X' 'FinishBroken.sh ended - Proxima-like system!'
	else
		./email.sh $1 'X/X' 'FinishBroken.sh ended'
	fi	

