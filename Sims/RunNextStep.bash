#!/bin/sh
############################################################################### 
### Run the next step of a simulation, with the stop time specified exactly
### e.g: ./RunNextStep.bash dir >> dir/run.pipe &

### Start clock
t1=$(date +%s)


### Five stop times = [2.8, 4.6, 6.4, 8.2, 10]*365.25e7
stop=2.8*365.25e7
echo $stop

pwd=$(pwd)

echo '================================================================'
echo 'Continuing simulation: '$1', t='$stop

### Update param.dmp file
python -c 'import AlphaCenModule as AC;AC.WriteParam("'$1'/Out/param.dmp", '$stop', mA=.123)'

### Copy most recent completed run to backup
if [ $k -gt $mintime ]; then
	\cp -p $1/Out/*.out $1/Out/Backup
	\cp -p $1/Out/*.dmp $1/Out/Backup
	\cp -p $1/Out/*.tmp $1/Out/Backup
	echo 'Attempted backup of previous timestep'
fi

### Run the next step
cd $1/Out
	./merc_disk
	./elem
	mv *.aei AeiOutFiles
cd $pwd

# For long simulations, write looptime
	# Stop clock for iteration
	t2=$(date +%s)
echo $dir'	'$stop'	'$(echo "$t2 - $t1"|bc ) >> disktimes.txt


### Run summary (pick appropriate masses)
python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", '$stop'/365.25, WhichTime="****", wantsum=False, mA=.123, mB=.123, mC=.123)'

R CMD BATCH -$dir '../Analysis/ReadDisk.R'

