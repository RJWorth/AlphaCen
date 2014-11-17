#!/bin/sh
############################################################################### 
### Run the next step of a simulation, with the stop time specified exactly
### e.g: ./RunNextStep.bash dir >> dir/run.pipe &

### Five stop times = [2.8, 4.6, 6.4, 8.2, 10]*365.25e7
stop=1.01*365.25e7
echo $stop

mA=0.123

### Start clock
t1=$(date +%s)

pwd=$(pwd)

echo '================================================================'
echo 'Continuing simulation: '$1', t='$stop

### Copy most recent completed run to backup
if [ $k -gt $mintime ]; then
	\cp -p $1/Out/*.out $1/Out/Backup
	\cp -p $1/Out/*.dmp $1/Out/Backup
	\cp -p $1/Out/*.tmp $1/Out/Backup
	echo 'Attempted backup of previous timestep'
fi

### Update param.dmp file
python -c 'import AlphaCenModule as AC;AC.WriteParam("'$1'/Out/param.dmp", '$stop', mA=.123)'

### Run the next step
cd $1/Out
	./merc_disk
	./elem
	mv *.aei AeiOutFiles
cd $pwd

# For long simulations, write looptime
	# Stop clock for iteration
	t2=$(date +%s)
echo $1'	'$stop'	'$(echo "$t2 - $t1"|bc ) >> disktimes.txt


### Run summary (pick appropriate masses)
if [ ${1: -1:1} = 2 ]; then
	python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", '$stop', WhichTime="Disk", wantsum=False, wantplot=False, mode="binary", mA='$mA', mB=.123)'
else
	python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", '$stop', WhichTime="Disk", wantsum=False, wantplot=False, mode="triple", mA='$mA', mB=.123)'
fi

R CMD BATCH -$1 '../Analysis/ReadDisk.R'

