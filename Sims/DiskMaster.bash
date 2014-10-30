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
Which=()
dir='/Disk1'

### Simulation parameters
mA=0.123	#1.105
mB=0.123	#0.934
mC=0.123

newdisk=T	# generate a new disk, T or F
amin=0.1	# minimum extent of disk, in AU
sz=small   	# add to big.in or small.in? (not used yet)

############### Set up simulations
for i in ${Which[*]}
do
	j=$sim$i$dir

	# Set up triple system sim
	if [ ! -d $j-C ]; then
		echo 'Creating '$j-C
		\cp -rp $sim$i/Original $j-C
	else
		echo $j-C exists
	fi
	# If no backup dir yet, make one
	mkdir -p $j-C/Out/Backup
	# Reset stop files
	echo 'False' > $j-C/stopfile.txt
	echo 'False' > $j-C/bigstopfile.txt

	# Generate disk
	if [ $newdisk = T ]; then
		python -c 'import AlphaCenModule as AC; AC.MakeSmallTestDisk("'$j'-C",amin='$amin',objs=[],size="'$sz'")'
	fi
	
	# Set up binary system sim
#	if [ ! -d $j-B ]; then
#		echo 'Creating '$j-B
#		\cp -rp $j-C $j-B
#	else
#		\cp -p $j-C/In/small.in $j-B/In/small.in
#		mkdir -p $j-C/Out/Backup
#		echo 'False' > $j-C/stopfile.txt
#		echo 'False' > $j-C/bigstopfile.txt
#		echo $j-B exists, copying small list over
#	fi
	if [ -d $j-B ]; then
		\rm -r $j-B
	fi
	echo 'Creating '$j-B
	\cp -rp $j-C $j-B

	head -10 $j-C/In/big.in > $j-B/In/big.in

	# Set up single star sim
#	if [ ! -d $j-A ]; then
#		echo 'Creating '$j-A
#		\cp -rp $j-C $j-A
#	else
#		\cp -p $j-C/In/small.in $j-A/In/small.in
#		echo $j-A exists, copying small list over
#	fi
#	head -6 $j-C/In/big.in > $j-A/In/big.in

	# Start all sims
	nice -n 10 ./DiskRun.bash $j-C $mA T > $j-C/run.pipe &
	echo 'master: '$j'-C  '$!
	nice -n 10 ./DiskRun.bash $j-B $mA T > $j-B/run.pipe &
	echo 'master: '$j'-B  '$!
#	nice -n 10 ./DiskRun.bash $j-A $mA T > $j-A/run.pipe &
#	echo 'master: '$j'-A  '$!

done

