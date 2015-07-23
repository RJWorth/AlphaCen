#!/bin/bash
##############################################################################
# This script cleanly reruns the most recent Merc95
t1=$(date +%s)

machine=$(hostname -s)
dir=$1
ti=$2	# = log(years)
tf=$3	# = log(years)

c=../Code
in=In/
out=Out/

mStar=0.934
rEj=1e5

### Range for iterations
if [ $machine = chloe ]; then
	timerange=$(jot $(echo "$tf-$tie+1" | bc) $mintime)
else
	timerange=$(seq $ti $tf)
fi
echo '	timerange '$timerange

### Move to working dir
cd $dir

# write param.in file
python -c "import Merc as M; M.WriteParamInFile(loc='"$in"', tf=365.25e"$ti", mStar="$mStar", rEj="$rEj")"

# write big.in file with stars+planetesimals (currently only has planetesimals around AlCenA)
python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(), loc='"$in"')"

### Compile and run
gfortran -w -o merc_$1 $c/mercury6_2.f95 $c/drift.f95 $c/orbel.f95 $c/mal.f95 $c/mce.f95 $c/mco.f95 $c/mdt.f95 $c/mio.f95 $c/mfo.f95 $c/mxx.f95 $c/both.f95
gfortran -w -o elem_$1 $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95

### Clean previous sim
rm Out/*
rm Aei/*

echo 'run: '$1', t = '$timerange
echo '==================== start t = 1e'$ti'-'$tf' yrs ===================='
	### Loop over time lengths
		for k in $timerange; do
		# Write param.dmp file
		if [ $k != $ti ]; then
		python -c "import Merc as M; M.WriteParamInFile(loc='"$out"', f='dmp', tf=365.25e"$k", mStar="$mStar", rEj="$rEj")"
		fi
		#### Run mercury
#		cd $1/Out;	./merc_$end;	cd $pwd
		./merc_$1
		echo '----------------------- t = 1e'$k' yrs -------------------------'
		done
echo '====================== end t = 1e'$ti'-'$tf' yrs ===================='


### Run element
./elem_$1

### Organize files
mv *.tmp Out
mv *.aei Aei

### Print runtime
t2=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t2'))'

cd ..

dt=$(python -c 'print('$t2'-'$t1')')
echo $1'    '$2'    '$dt >> runtimes.txt


