#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir

# Which sim(s)?
newrun=T

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

### Set up triple system sim

	if [ ! -d $j'A-3' ]; then
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
	
	# Start triple sim
	nice -n 10 ./DiskRun.bash $j'A-3' $mA $newrun > $j'A-3'/run.pipe &
	echo 'master: '$j'A-3  '$!

### Set up binary system sim

	if [ -d $j'A-2' ]; then
		\rm -r $j'A-2'
	fi
	echo 'Creating '$j'A-2'
	\cp -rp $j'A-3' $j'A-2'

	head -10 $j'A-3'/In/big.in > $j'A-2'/In/big.in


	# Start binary sim
	nice -n 10 ./DiskRun.bash $j'A-2' $mA $newrun > $j'A-2'/run.pipe &
	echo 'master: '$j'A-2  '$!

done

