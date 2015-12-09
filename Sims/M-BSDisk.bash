#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir
machine=$(hostname -s)
home=$(pwd)

# Which sim(s)?
newrun=T

#sim='Proxlike/Prx'
sim='BS/Prx'
dir='/Disk'
### Assign directories based on machine
#if [ ${machine:0:5} = 'lionx' ]; then
#	echo 'running on '$machine
#	if [ ${machine:5:1} = 'j' ]; then
#	Which=(01 02 03 04 05 06 07)
#	elif  [ ${machine:5:1} = 'g' ]; then
#	Which=(08 09 10 11 12 13 14)
#	elif  [ ${machine:5:1} = 'f' ]; then
#	Which=(15 16 17 18 19 20 21)
#	elif  [ ${machine:5:1} = 'i' ]; then
#	Which=(22 23 24 25 26 27 28)
#	elif  [ ${machine:5:1} = 'h' ]; then
#	Which=(29 30 31 32 33 34)
#	else
#		echo 'unknown lionx machine'
#	fi
#else
#	echo 'Non-lionx machine'
#	Which=()
#fi
#Which=(01 03 05 06 07 09 10 11 12 13 14 15 16 17 19 20 21 22 23 24 26 33) #
#Which=(02 04 08 18 25 27 28 29 30 31 32 34)
#Which=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)
Which=(34)
echo ${Which[*]}

### Simulation parameters
mA=0.123	#1.105
mB=0.123	#0.934
mC=0.123

newdisk=F	# generate a new disk, T or F
#amin=0.1	# minimum extent of disk, in AU
#sz=small   	# add to big.in or small.in? (not used yet)

############### Set up simulations
for i in ${Which[*]}
do
	j=$sim$i$dir

### Set up triple system sim (disk around A)
#	cent='AlCenA'
#	if [ $newrun = T ]; then
#		if [ ! -d $j'A-3' ]; then
#			echo 'Creating '$j'A-3'
#			\cp -rp $sim$i/Original $j'A-3'
#		else
#			echo $j'A-3' exists
#		fi
 #               \rm $j'A-3'/run.pipe
#		# If no backup dir yet, make one
#		mkdir -p $j'A-3'/Out/Backup
#		# Reset stop files
#		echo 'False' > $j'A-3'/stopfile.txt
#		echo 'False' > $j'A-3'/bigstopfile.txt
#
#		# Generate disk
#		if [ $newdisk = T ]; then
#		python -c 'import AlphaCenModule as AC; AC.MakeSmallTestDisk("'$j'A-3",centobj="'$cent'",m=0.0)'
#		fi
#	fi	

### Set up binary system sim (disk around A)
#	if [ $newrun = T ]; then
#		if [ -d $j'A-2' ]; then
#			\rm -r $j'A-2'
#		fi
#		echo 'Creating '$j'A-2'
#		\cp -rp $j'A-3' $j'A-2'
#
#		head -10 $j'A-3'/In/big.in > $j'A-2'/In/big.in
#	fi

### Move data into 'Original' dir
	if [ ! -d $sim$i/Original ]; then
		echo 'No Original directory; creating one now'
		mkdir $sim$i/Original
		mv $sim$i/* $sim$i/Original
	fi

### Set up triple system sim (disk around B)
	cent='AlCenB'
	if [ $newrun = T ]; then
		if [ ! -d $j'B-3' ]; then
			echo 'Creating '$j'B-3'
			\cp -rp $sim$i/Original $j'B-3'
		else
			echo $j'B-3' exists
		fi
		\rm ${j}B-3/run.pipe
		# If no backup dir yet, make one
		mkdir -p $j'B-3'/Out/Backup
		# Reset stop files
		echo 'False' > $j'B-3'/stopfile.txt
		echo 'False' > $j'B-3'/bigstopfile.txt

		# Generate disk
		if [ $newdisk = T ]; then
		python -c 'import AlphaCenModule as AC; AC.MakeSmallTestDisk("'$j'B-3",centobj="'$cent'",m=0.0)'
		fi
	else
		mkdir -p $j'B-3'/Out/Backup
		\cp -p $j'B-3'/Out/*.out    $j'B-3'/Out/Backup
		\cp -p $j'B-3'/Out/*.dmp    $j'B-3'/Out/Backup
		\cp -p $j'B-3'/Out/*.tmp    $j'B-3'/Out/Backup

#		\cp -p $j'B-3'/Out/Backup/* $j'B-3'/Out
		echo 'Attempted backup of previous timestep'
	fi

### Set up binary system sim (disk around B)
	if [ $newrun = T ]; then
		if [ -d $j'B-2' ]; then
			\rm -r $j'B-2'
		fi
		echo 'Creating '$j'B-2'
		\cp -rp $j'B-3' $j'B-2'

		head -10 $j'B-3'/In/big.in > $j'B-2'/In/big.in
	else
		mkdir -p $j'B-2'/Out/Backup
		\cp -p $j'B-2'/Out/*.out    $j'B-2'/Out/Backup
		\cp -p $j'B-2'/Out/*.dmp    $j'B-2'/Out/Backup
		\cp -p $j'B-2'/Out/*.tmp    $j'B-2'/Out/Backup

#		\cp -p $j'B-2'/Out/Backup/* $j'B-2'/Out
		echo 'Attempted backup of previous timestep'
	fi

	# Start running triple and binary sims
	if [ $machine = 'shapiro' ] || [ $machine = 'chloe' ] || [ ${machine:0:6} = 'hammer' ]; then
		echo 'using bash script'
#		nice -n 10 ./R-BSDisk.bash $j'A-2' $mA $newrun > $j'A-2'/run.pipe &
#		echo 'master: '$j'A-2  '$!
#		nice -n 10 ./R-BSDisk.bash $j'A-3' $mA $newrun > $j'A-3'/run.pipe &
#		echo 'master: '$j'A-3  '$!
		nice -n 10 ./R-BSDisk.bash $j'B-2' $mA $newrun > $j'B-2'/run.pipe &
		echo 'master: '$j'B-2  '$!
		nice -n 10 ./R-BSDisk.bash $j'B-3' $mA $newrun > $j'B-3'/run.pipe &
		echo 'master: '$j'B-3  '$!
	elif [ ${machine:0:5} = 'lionx' ]; then
		echo 'using qsub script'
#		d=${j}A-2
#	        runname=${d:12:2}_${d:19:3}
#		qsub -v h=$home,dir=$d,mA=$mA,newrun=$newrun,runname=$runname run1.pbs
#               d=${j}A-3
#               runname=${d:12:2}_${d:19:3}
#		qsub -v h=$home,dir=$d,mA=$mA,newrun=$newrun,runname=$runname run1.pbs
                d=${j}B-2
                runname=${d:12:2}_${d:19:3}
		if [ $newrun = T ]; then
			qsub -v h=$home,dir=$d,mA=$mA,newrun=$newrun run1.pbs
		elif [ $newrun = F ]; then
			qsub -v h=$home,dir=$d,mA=$mA,mB=$mB,newrun=$newrun run2.pbs
		fi
                d=${j}B-3
                runname=${d:12:2}_${d:19:3}
		if [ $newrun = T ]; then
			qsub -v h=$home,dir=$d,mA=$mA,newrun=$newrun run1.pbs
		elif [ $newrun = F ]; then
			qsub -v h=$home,dir=$d,mA=$mA,mB=$mB,newrun=$newrun run2.pbs
		else
			echo 'invalid newrun parameter, not running!'
		fi 
	else
		echo 'unknown host machine -- not running!!!'
	fi
done

