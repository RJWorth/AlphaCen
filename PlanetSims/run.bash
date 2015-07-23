#!/bin/bash
##############################################################################
# This script cleanly reruns the most recent Merc95
t1=$(date +%s)

c=../Code
in=In/

cd $1

### Clean previous sim
rm Out/*
rm Aei/*

# write param.in file
python -c "import Merc as M; M.WriteParamInFile(loc='"$in"', tf=365.25e"$2", mStar=1.105, rEj=100, rStar=0.005)"

# write big.in file with stars+planetesimals (currently only has planetesimals around AlCenA)
python -c "import Merc as M; M.WriteObjInFile(objlist=M.GetObjList(), loc='"$in"')"

### Compile and run
gfortran -w -o merc $c/mercury6_2.f95 $c/drift.f95 $c/orbel.f95 $c/mal.f95 $c/mce.f95 $c/mco.f95 $c/mdt.f95 $c/mio.f95 $c/mfo.f95 $c/mxx.f95 $c/both.f95
./merc
gfortran -w -o elem $c/element6.f95 $c/e_sub.f95 $c/orbel.f95 $c/both.f95
./elem

### Organize files
mv *.tmp Out
mv *.aei Aei

### Print runtime
t2=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t2'))'

cd ..

dt=$(python -c 'print('$t2'-'$t1')')
echo $1'    '$2'    '$dt >> runtimes.txt


