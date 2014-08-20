#!/bin/sh
############################################################################### 
# Start time
t1=$(date +%s)
machine=$(hostname -s)

### Simulation parameters

cent=A		# placeholder; this currently does nothing

mA=0.123	# mSun, actual = 1.105
mB=0.123	# mSun, actual = 0.934
mC=0.123	# mSun, actual = 0.123

aBmin=23.4	# AU
aBmax=28.1
eBmin=0.0	# 0 to 1
eBmax=0.52
iBmin=0.0	# degrees
iBmax=0.0

aCmin=1.0	# *aB(1+eB)/(1-eC), i.e. min(periapse of C) = this*(apoapse of B)
aCmax=10.0	# ditto
eCmin=0.0
eCmax=0.75
iCmin=0.0
iCmax=180.0

vers='ury_TG.for'	# merc+vers=filename for mercury
mintime=3	# = log(years)
maxtime=9	# = log(years)
output=3	# = log(years)
step=10.0	# = days
niter=5 	# = number of iterations to run
user='yes'	# use user-defined forces?

### Write files.in
./writefiles.bash $1 $vers

### Range for iterations
if [ $machine = chloe ]; then
	itrange=$(jot $niter 1)	# start at y, go up x-1 steps?
	timerange=$(jot $(echo "$maxtime-$mintime+1" | bc) $mintime)
else
	itrange=$(seq 1 $niter)	# go from x to y
	timerange=$(seq $mintime $maxtime)
fi
echo '	itrange '$itrange
echo '	timerange '$timerange
### Do iterations
	for j in $itrange
	do
	echo 'run: '$1', '$j
	# Start clock for iteration
	t3=$(date +%s)
	# Clean out old sim
	\rm $1/Out/*.dmp
	\rm $1/Out/*.out
	#### Randomize B and C parameters (a, e, i)
	python -c 'import AlphaCenModule; AlphaCenModule.MakeBigRand( "'$1'",'$j', "'$cent'", 	'$aBmin','$aBmax','$eBmin','$eBmax','$iBmin','$iBmax', '$aCmin','$aCmax','$eCmin','$eCmax','$iCmin','$iCmax', mB='$mB',mC='$mC')'

	# Write param.in file
	./writeparam.bash $1 $mintime $output $step $mintime $user $mA
	# Compile mercury
	gfortran -w -O1 -o $1/Out/merc_AC$1 Files/merc$vers
		#j in, to fix colors
	### Loop over time lengths
	for k in $timerange; do

		#### Run mercury
		cd $1/Out;	./merc_AC$1;	cd ../..

		### Compile and run Element to get orbit details
		if [ $machine = chloe ]; then	
			gfortran-4.2 -w -O1 -o $1/Out/elem Files/elem.for 	
				#j in, to fix colors
		else
			gfortran -w -O1 -o $1/Out/elem Files/elem.for 		
				#j in, to fix colors
		fi
		cd $1/Out;	./elem;	cd ../..
		\mv $1/Out/*.aei $1/Out/AeiOutFiles
		### Summarize iteration; write if stop conditions reached
		python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'", 1e'$k', 1e'$maxtime', WhichTime="'$j'", cent="'$cent'", machine="'$machine'", wantsum=True, wantplot=False, mode="triple", mA='$mA', mB='$mB', mC='$mC')'

		# For long simulations, write looptime
		if [ $k -ge 7 ]; then
			# Get end time
			endtime=$(python -c 'import datetime;print(datetime.datetime.now())')
			# Stop clock for iteration
			t4=$(date +%s)
			echo $1'	'$machine'	'$j'	'$k'	'${endtime:0:16}'	'$(echo "$t4 - $t3"|bc ) >> looptimes.txt
		fi	# k>=7
		# If stopfile=true, don't continue this simulation
		stop=$(cat $1/stopfile.txt)
		if [ $stop = 'True' ]; then
			break
		else
		# Write param.dmp file for next timestep
		./writeparam.bash $1 $(echo "$k+1"|bc) $output $step $mintime $user $mA
		fi	# stop
	done	#timerange
	# Check for bigstop flag to stop the 'niter' loop
	bigstop=$(cat $1/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Proxima-like system?'
		./email.sh $1 $j'/'$niter 'Proxima-like system?'
		break
	fi
	# Check continuation flag to stop the 'niter' loop
	cont=$(cat ContinuationFlag.txt)
	if [ $cont = 'False' ]; then
		echo 'Iterations stopped by continuation flag.'
		./email.sh $1 $j'/'$niter 'Iterations stopped by continuation flag'
		break
	fi
done	#niter

# Write stop time for this directory:
t2=$(date +%s)

echo $1"	"$machine"	"$j"/"$niter"	"$user"	"$vers"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt

./email.sh $1 $j'/'$niter 'AC finished'

