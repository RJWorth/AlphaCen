#!/bin/bash
set -e	# stop script on error
##############################################################################
# This script cleanly reruns the most recent Merc95
t1=$(date +%s)

machine=$(hostname -s)
dir=$1
ti=$2	# = log(years)
tf=$3	# = log(years)
rtr=$4  # truncation radius (AU)
sigma=$5 # density scale factor

home=..
c=../Code
in=In/
out=Out/

mStar=0.934
rEj=1e2

### Range for iterations
if [ $machine = chloe ]; then
	timerange=$(jot $(echo "$tf-$ti+1" | bc) $ti)
else
	timerange=$(seq $ti $tf)
fi
echo '	timerange '$timerange

### Move to working dir
cd $dir

# write param.in file
tOut=$(python -c 'print('$ti'-3)')
python -c "import Merc as M; M.WriteParamInFile(loc='"$in"', tf=365.25e"$ti", tOut=365.25e"$tOut", mStar="$mStar", rEj="$rEj")"

# write big.in file with stars+planetesimals (currently only has planetesimals around AlCenA)
#python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(rtr="$rtr",sigC="$sigma",iMax=0.), loc='"$in"')"

### Compile and run
gfortran -w -o merc_$1 $c/mercury6_2.f95 $c/drift.f95 $c/orbel.f95 $c/mal.f95 $c/mce.f95 $c/mco.f95 $c/mdt.f95 $c/mio.f95 $c/mfo.f95 $c/mxx.f95 $c/both.f95
if [ $machine = chloe ]; then
	gfortran-4.2 -w -o elem_$1 $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95
else
	gfortran -w -o elem_$1 $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95
fi

### Clean previous sim
if [ $(ls -l Out | wc -l) -gt 0 ]; then rm Out/*; fi
if [ $(ls -l Aei | wc -l) -gt 0 ]; then rm Aei/*; fi

echo 'run: '$1', t = '$timerange
echo '==================== start t = 1e'$ti'-'$tf' yrs ===================='
	### Loop over time lengths
		for k in $timerange; do

	### Write param.dmp file
		if [ $k != $ti ]; then
		tOut=$(python -c 'print('$k'-3)')
		python -c "import Merc as M; M.WriteParamInFile(loc='"$out"', f='dmp', tf=365.25e"$k", tOut=365.25e"$tOut", mStar="$mStar", rEj="$rEj")"
		fi

	#### Run mercury
		./merc_$1
		echo '----------------------- t = 1e'$k' yrs -------------------------'

	### Backup
	#	\cp -p Out/* Backup

	### Print runtime thus far
		t2=$(date +%s)
		dt=$(python -c 'print('$t2'-'$t1')')
		echo $dir'    '$k'    '$dt >> $home/runtimes.txt

		done
echo '====================== end t = 1e'$ti'-'$tf' yrs ===================='


### Run element
./elem_$1

### Organize files
if [ $(ls -l *.tmp | wc -l) -gt 0 ]; then mv *.tmp Out; fi
if [ $(ls -l *.aei | wc -l) -gt 0 ]; then mv *.aei Aei; fi

### Print runtime
t3=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t3'))'

cd $home

dt=$(python -c 'print('$t3'-'$t1')')
echo $dir'    '$k'    '$dt >> runtimes.txt


