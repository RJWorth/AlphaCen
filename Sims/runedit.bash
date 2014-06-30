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

vers='ury_TG_edit.for'	# merc+vers=filename for mercury
maxtime=10	# = log(years)
mintime=10	# = log(years)
output=3	# = log(years)
step=10.0	# = days
niter=1		# = number of iterations to run
user='yes'	# use user-defined forces?

### Write files.in
#./writefiles.sh $1		# for .for versions
./writefiles.bash $1	# for debugging.f vers

### Range for iterations
if [ $machine = chloe ]; then
	itrange=$(jot $niter 1)	# start at y, go up x-1 steps?
	timerange=$(jot $(echo "$maxtime-$mintime+1" | bc) $mintime)
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
	#	python -c 'import AlphaCenModule; AlphaCenModule.MakeBigRand( "'$1'",'$j',"'$cent'", 	'$aBmin','$aBmax','$eBmin','$eBmax','$iBmin','$iBmax', '$aCmin','$aCmax','$eCmin','$eCmax','$iCmin','$iCmax')'

	# Write param.in file
	./writeparam.bash $1 $mintime $output $step $mintime $user
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
		#python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'","'$j'","'$cent'",1e'$maxtime',1e'$k')'

		# If stopfile=true, don't continue this simulation
		stop=$(cat $1/stopfile.txt)
		#if [ $stop = 'True' ]; then
		#	break
		#fi

		# Write param.dmp file for next timestep
		./writeparam.bash $1 $(echo "$k+1"|bc) $output $step $mintime $user
	done	#timerange
	#bigstop=$(cat $1/bigstopfile.txt)
	#if [ $bigstop = 'True' ]; then
	#	echo 'Bigstop reached! Prox-like system.'
	#	#./email.bash $1" prox-like"
	#	break
	#fi
done	#niter

# Write stop time for this directory:
t2=$(date +%s)
echo $1"	"$machine"	"$niter"	"$user"	"$(echo "$t2 - $t1"|bc )"	"$vers >> testtimes.txt
#./email.bash $1" done"

