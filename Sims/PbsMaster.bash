#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir

# Which sim(s)?
newrun=F

sim='Proxlike/Prx'
Which=(01)
#Which=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
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
	if [ $newrun = T ]; then
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
	fi	
	# Start running triple system
	qsub -v dir=$j-C,mA=$mA,newrun=$newrun -o $j-C/run.pipe -j oe run1.pbs

	# Set up binary system sim
	if [ $newrun = T ]; then
		if [ -d $j-B ]; then
			\rm -r $j-B
		fi
		echo 'Creating '$j-B
		\cp -rp $j-C $j-B

		head -10 $j-C/In/big.in > $j-B/In/big.in
	fi

	# Start binary system
	qsub -v dir=$j-B,mA=$mA,newrun=$newrun -o $j-B/run.pipe -j oe run1.pbs

done

