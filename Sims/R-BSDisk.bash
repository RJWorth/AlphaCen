#!/bin/bash
###############################################################################
# Run the simulation in directory $1, with parameters set below
# > ./R-BSDisk.bash BS/Prx#/Disk# BSDisk.inc > Dir/run.pipe &

# include file
inc=$2
source $inc

# Start time
t1=$(date +%s)
machine=$(hostname -s)
home=$(pwd)
echo $home

### Simulation parameters
dir=$1
tag=${dir:6:2}

# Compile mercury and element
#if [ $newrun = T ]; then
#	if [ $newmerc = T ]; then
		gfortran -w -O1 -o $dir/Out/merc_$tag Files/merc$vers
#	fi
#fi
	gfortran -w -O1 -o $dir/Out/elem_$tag Files/elem.for 		
		#j in, to fix colors
	
### Do iterations
	echo 'DiskRun: '$dir', '$j
	# Start clock for iteration
	t3=$(date +%s)
	echo '================================================================'
	# Clean out old sim
	if [ $newrun = T ]; then
		\rm $dir/Out/*.dmp
		\rm $dir/Out/*.out
		echo 'cleaning old data'
	else
		echo 'continuing previous run'
	fi

	# Write param.in file
	if [ $newrun = T ]; then
		./writeparam.bash $dir $mintime $output $step $mintime $user $mA $alg
	else
		./writeparam.bash $dir $mintime $output $step        0 $user $mA $alg
	fi
	### Loop over time lengths
	for k in $timerange; do

		### Copy most recent completed run to backup
		if [ $newrun = T ]; then
			if [ $k -gt $mintime ]; then
			\cp -p $dir/Out/*.out $dir/Out/Backup
			\cp -p $dir/Out/*.dmp $dir/Out/Backup
			\cp -p $dir/Out/*.tmp $dir/Out/Backup
			echo 'Attempted backup of 1e'$(echo "$k-1"|bc )' timestep'
			fi
		fi

		#### Run mercury
		cd $dir/Out;	./merc_$tag;	cd $home

		### Run Element to get orbit details
		cd $dir/Out;	./elem_$tag;	cd $home
		\mv $dir/Out/*.aei $dir/Out/AeiOutFiles
		### Summarize iteration; write if stop conditions reached
		if [ ${dir: -1:1} = 2 ]; then
		python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$dir'", 1e'$k', WhichTime="Disk", wantsum=False, wantplot=False, mode="binary", mA='$mA', mB=.123, mC=.123)'
		else
		python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$dir'", 1e'$k', WhichTime="Disk", wantsum=False, wantplot=False, mode="triple", mA='$mA', mB=.123, mC=.123)'
		fi

#		if [ $k -ge 5 ]; then
#			R CMD BATCH -$dir '../Analysis/ReadDisk.R'
#		fi	# k >= 5

		# For long simulations, write looptime
		if [ $k -ge 3 ]; then
			# Get end time
			endtime=$(python -c 'import datetime;print(datetime.datetime.now())')
			size=$(python -c 'import os;print(os.path.getsize("'$dir'/Out/xv.out"))')
			# Stop clock for iteration
			t4=$(date +%s)
			echo $dir'	'$k'	'$(echo "$t4 - $t3"|bc ) >> disktimes.txt
		fi	# k>=7
		# If stopfile=true, don't continue this simulation
		stop=$(cat $dir/stopfile.txt)
		if [ $stop = 'True' ]; then
			break
		elif [ $k -lt $maxtime ]; then
		# Write param.dmp file for next timestep
		echo '----------------------------------------------------------------'
		./writeparam.bash $dir $(echo "$k+1"|bc) $output $step $mintime $user $mA $alg
		fi	# stop
	done	#timerange
	# Check for bigstop flag to stop the 'niter' loop
	bigstop=$(cat $dir/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Proxima-like system?'
#		./email.sh $dir $j'/'$niter 'Disk sim completed'
		break
	fi
	# Check continuation flag to stop the 'niter' loop
	cont=$(cat ContinuationFlag.txt)
	if [ $cont = 'False' ]; then
		echo 'Iterations stopped by continuation flag.'
#		./email.sh $dir $j'/'$niter 'Iterations stopped by continuation flag'
		break
	fi

#R CMD BATCH -$dir '../Analysis/ReadDisk.R'

# Write stop time for this directory:
t2=$(date +%s)


#echo $dir"	"$machine"	"$j"/"$niter"	"$user"	"$vers"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt

#./email.sh $dir $j'/'$niter 'Prx disk run finished'

