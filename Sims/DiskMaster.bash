#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir

# to clean previous run:
# move plots, SumAll, AllParams, Rout files to Saved/whatever
# \cp -p BlankDir/summary.out *Dir*		or
# \rm *Dir*/summary.out
# \rm *Dir*/InParams.txt

# Which sim(s)?
sim='Proxlike/Prx'
Which=(01)
dir='/Disk1'

### Simulation parameters
mA=0.123	#1.105
mB=0.123	#0.934
mC=0.123

newdisk=T	# generate a new disk, T or F
amin=0.1	# minimum extent of disk, in AU
sz=small   	# add to big.in or small.in?

############### Set up simulations
for i in ${Which[*]}
do
	j=$sim$Which$dir	

	# Set up triple system sim
	if [ ! -d $j-C ]; then
		echo 'Creating '$j-C
		\cp -rp $sim$Which/Original $j-C
	else
		echo $j-C exists
	fi

	# Generate disk
	if [ $newdisk = T ]; then
		python -c 'import AlphaCenModule as AC; AC.MakeSmallTestDisk("'$j'-C",amin='$amin',objs=[],size="'$sz'")'
	fi
	
	# Set up binary system sim
	if [ ! -d $j-B ]; then
		echo 'Creating '$j-B
		\cp -rp $j-C $j-B
	else
		\cp -p $j-C/In/small.in $j-B/In/small.in
		echo $j-B exists, copying small list over
	fi
	head -10 $j-C/In/big.in > $j-B/In/big.in

	# Set up single star sim
	if [ ! -d $j-A ]; then
		echo 'Creating '$j-A
		\cp -rp $j-C $j-A
	else
		\cp -p $j-C/In/small.in $j-A/In/small.in
		echo $j-A exists, copying small list over
	fi
	head -6 $j-C/In/big.in > $j-A/In/big.in

	# Start all sims
	nice -n 10 ./DiskRun.bash $j-C > $j-C/run.pipe &
	echo 'master: '$j'-C  '$!
	nice -n 10 ./DiskRun.bash $j-B > $j-B/run.pipe &
	echo 'master: '$j'-B  '$!
	nice -n 10 ./DiskRun.bash $j-A > $j-A/run.pipe &
	echo 'master: '$j'-A  '$!

done

