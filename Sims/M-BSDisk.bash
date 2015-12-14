#!/bin/bash
############################################################################### 
# Start an instance of run.sh in each Dir
machine=$(hostname -s)
home=$(pwd)

### Include filename
inc=BSDisk.inc
source $inc

############### Set up simulations
for i in ${DirNums[*]}
do
	j=$dir1$i$dir2

### Set up triple system sim (disk around A)
### Set up binary system sim (disk around A)

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
		nice -n 10 ./R-BSDisk.bash $j'B-2' $inc > $j'B-2'/run.pipe &
		echo 'master: '$j'B-2  '$!
		nice -n 10 ./R-BSDisk.bash $j'B-3' $inc > $j'B-3'/run.pipe &
		echo 'master: '$j'B-3  '$!
	elif [ ${machine:0:5} = 'disable' ]; then
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

