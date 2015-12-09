#!/bin/bash
set -e	# stop script on error
##############################################################################
# This script cleanly reruns the most recent Merc95
t1=$(date +%s)

dir=$1
source $2

### Range for iterations
if [ $machine = chloe ]; then
	timerange=$(jot $(echo "$tf-$ti+1" | bc) $ti)
else
	timerange=$(seq $ti $tf)
fi
echo '	timerange '$timerange

### Move to working dir
cd $dir

### Write param.in file
tOut=$(python -c 'print('$k'-'$tDmp')')
python -c "import Merc as M; M.WriteParamInFile(loc='"$in"', tf=365.25e"$ti", tOut=365.25e"$tOut", mStar="$mStar", rEj="$rEj")"

### Write big.in file with stars+planetesimals (currently only has planetesimals around AlCenA)
if [ $new = T ]; then
	python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(rtr="$rtr",sigC="$sigma",starA='$starA',iMax='$iMax'), loc='"$in"')"
fi

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
		tOut=$(python -c 'print('$k'-'$tDmp')')
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


