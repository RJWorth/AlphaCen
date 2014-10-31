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
dir='/Disk'

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
		echo 'Creating '$j'A-3'
		\cp -rp $sim$i/Original $j'A-3'
	else
		echo $j'A-3' exists
	fi
	# If no backup dir yet, make one
	mkdir -p $j'A-3'/Out/Backup
	# Reset stop files
	echo 'False' > $j'A-3'/stopfile.txt
	echo 'False' > $j'A-3'/bigstopfile.txt

	# Generate disk
	if [ $newdisk = T ]; then
		python -c 'import AlphaCenModule as AC; AC.MakeSmallTestDisk("'$j'A-3",amin='$amin',objs=[],size="'$sz'")'
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
	if [ -d $j'A-2' ]; then
		\rm -r $j'A-2'
	fi
	echo 'Creating '$j'A-2'
	\cp -rp $j'A-3' $j'A-2'

	head -10 $j'A-3'/In/big.in > $j'A-2'/In/big.in

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
	nice -n 10 ./DiskRun.bash $j'A-3' $mA T > $j'A-3'/run.pipe &
	echo 'master: '$j'A-3  '$!
	nice -n 10 ./DiskRun.bash $j'A-2' $mA T > $j'A-2'/run.pipe &
	echo 'master: '$j'A-2  '$!
#	nice -n 10 ./DiskRun.bash $j-A $mA T > $j-A/run.pipe &
#	echo 'master: '$j'-A  '$!

done

