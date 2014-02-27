#!/bin/sh
############################################################################### 
# Start time
t1=$(date +%s)
machine=$(hostname -s)

### Simulation parameters

cent=A

aBmin=23.4
aBmax=28.1
eBmin=0.0
eBmax=0.52
iBmin=0.0
iBmax=0.0

aCmin=1.0	# *aB(1+eB)/(1-eC), i.e. min(periapse of C) = this*(apoapse of B)
aCmax=10.0	# ditto
eCmin=0.0
eCmax=0.75
iCmin=0.0
iCmax=180.0

maxtime=10	# = log(years)
mintime=5	# = log(years)
output=3	# = log(years)
step=10.0	# = days
niter=10	# = number of iterations to run
user='yes'	# use user-defined forces?

### Write files.in
./writefiles.sh $1

### Range for iterations
if [ $machine = chloe ]; then
	itrange=$(jot $niter 1)	# start at y, go up x-1 steps?
	timerange=$(jot $maxtime $mintime)
else
	itrange=$(seq 1 $niter)	# go from x to y
	timerange=$(seq $mintime $maxtime)
fi

### Do iterations
	for j in $itrange
	do
	echo 'run: '$1', '$j
	\rm $1/Out/*.dmp
	\rm $1/Out/*.out
	#### Randomize B and C parameters (a, e, i)
		python -c 'import AlphaCenModule; AlphaCenModule.MakeBigRand( "'$1'",'$j',"'$cent'", 	'$aBmin','$aBmax','$eBmin','$eBmax','$iBmin','$iBmax', '$aCmin','$aCmax','$eCmin','$eCmax','$iCmin','$iCmax')'

	# Write param.in file
	./writeparam.sh $1 $mintime $output $step $mintime $user
	# Compile mercury
	gfortran -w -o $1/Out/merc_AC$1 Files/mercury_TidesGas.for 
		#j in, to fix colors

	### Loop over time lengths
	for k in $timerange; do

		#### Run mercury
		cd $1/Out;	./merc_AC$1;	cd ../..

		### Compile and run Element to get orbit details
		if [ $machine = chloe ]; then	
			gfortran-4.2 -w -o $1/Out/elem Files/elem.for 	
				#j in, to fix colors
		else
			gfortran -w -o $1/Out/elem Files/elem.for 		
				#j in, to fix colors
		fi
		cd $1/Out;	./elem;	cd ../..
		\mv $1/Out/*.aei $1/Out/AeiOutFiles

		### Summarize iteration; write if stop conditions reached
		python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'","'$j'","'$cent'",1e'$maxtime',1e'$k')'

		# If stopfile=true, don't continue this simulation
		stop=$(cat $1/stopfile.txt)
		if [ $stop = 'True' ]; then
			break
		fi

		# Write param.dmp file for next timestep
		./writeparam.sh $1 $(echo "$k+1"|bc) $output $step $mintime $user
	done	#timerange
	bigstop=$(cat $1/bigstopfile.txt)
	if [ $bigstop = 'True' ]; then
		echo 'Bigstop reached! Prox-like system.'
		#./email.bash $1" prox-like"
		break
	fi
done	#niter

# Write stop time for this directory:
t2=$(date +%s)
echo $1"	"$machine"	"$niter"	"$user"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt
#./email.bash $1" done"

