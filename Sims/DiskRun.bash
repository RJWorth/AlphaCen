#!/bin/sh
###############################################################################
# Run the simulation in directory $1, with parameters set below
# > ./DiskRun.bash Proxlike/Prx#/Disk# > Dir/run.pipe &
# Start time
t1=$(date +%s)
machine=$(hostname -s)
home=$(pwd)
echo $home

### Simulation parameters
newmerc=T	# recompile the merc/elem executables?
vers='ury_TG.for'	# 'merc'+vers = filename for mercury

mintime=3	# = log(years)
maxtime=8	# = log(years)
output=1	# = log(years)
step=10.0	# = days
user='no'	# use user-defined forces?

mA=.123

### Range for iterations
if [ $machine = chloe ]; then
	timerange=$(jot $(echo "$maxtime-$mintime+1" | bc) $mintime)
else
	timerange=$(seq $mintime $maxtime)
fi
echo '	timerange '$timerange

# Compile mercury and element
if [ $newmerc = T ]; then
	gfortran -w -O1 -o $1/Out/merc_disk Files/merc$vers
	if [ $machine = chloe ]; then	
		gfortran-4.2 -w -O1 -o $1/Out/elem Files/elem.for 	
			#j in, to fix colors
	else
		gfortran -w -O1 -o $1/Out/elem Files/elem.for 		
			#j in, to fix colors
	fi
fi
### Do iterations
	echo 'DiskRun: '$1', '$j
	# Start clock for iteration
	t3=$(date +%s)
	echo '================================================================'
	# Clean out old sim
	\rm $1/Out/*.dmp
	\rm $1/Out/*.out

	# Write param.in file
#	./writeparam.bash $1 $mintime $output $step $mintime $user $mA
	### Loop over time lengths
	for k in $timerange; do

		#### Run mercury
		cd $1/Out;	./merc_disk;	cd $home

		### Run Element to get orbit details
		cd $1/Out;	./elem;	cd $home
		\mv $1/Out/*.aei $1/Out/AeiOutFiles
		### Summarize iteration; write if stop conditions reached
#		python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", 1e'$k', 1e'$maxtime', WhichTime="Disk", machine="'$machine'", wantsum=True, wantplot=False, mode="triple", mA='$mA', mB='$mB', mC='$mC')'

		if [ $k -ge 5 ]; then
			R CMD BATCH -$1 '../Analysis/ReadDisk.R'
		fi	# k >= 5

		# For long simulations, write looptime
		if [ $k -ge 1 ]; then
			# Get end time
			endtime=$(python -c 'import datetime;print(datetime.datetime.now())')
			size=$(python -c 'import os;print(os.path.getsize("'$1'/Out/xv.out"))')
			# Stop clock for iteration
			t4=$(date +%s)
			echo $1'	'$k'	'$(echo "$t4 - $t3"|bc ) >> disktimes.txt
		fi	# k>=7
		# If stopfile=true, don't continue this simulation
		stop=$(cat $1/stopfile.txt)
		if [ $stop = 'True' ]; then
			break
		else
		# Write param.dmp file for next timestep
		echo '----------------------------------------------------------------'
		./writeparam.bash $1 $(echo "$k+1"|bc) $output $step $mintime $user $mA
		fi	# stop
	done	#timerange
	# Check for bigstop flag to stop the 'niter' loop
	bigstop=$(cat $1/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Proxima-like system?'
#		./email.sh $1 $j'/'$niter 'Disk sim completed'
		break
	fi
	# Check continuation flag to stop the 'niter' loop
	cont=$(cat ContinuationFlag.txt)
	if [ $cont = 'False' ]; then
		echo 'Iterations stopped by continuation flag.'
		./email.sh $1 $j'/'$niter 'Iterations stopped by continuation flag'
		break
	fi

R CMD BATCH -$1 '../Analysis/ReadDisk.R'

# Write stop time for this directory:
t2=$(date +%s)


#echo $1"	"$machine"	"$j"/"$niter"	"$user"	"$vers"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt

./email.sh $1 $j'/'$niter 'Prx disk run finished'

