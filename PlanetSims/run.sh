#!/bin/bash
##############################################################################
# This script cleanly reruns the most recent Merc95
t1=$(date +%s)

l=../Code

cd $1

### Clean previous sim
rm Out/*
rm Aei/*

### Compile and run
gfortran -w -o merc $l/mercury6_2.f95 $l/drift.f95 $l/orbel.f95 $l/mal.f95 $l/mce.f95 $l/mco.f95 $l/mdt.f95 $l/mio.f95 $l/mfo.f95 $l/mxx.f95 $l/both.f95
./merc
gfortran -w -o elem $l/element6.f95 $l/e_sub.f95 $l/orbel.f95 $l/both.f95
./elem

### Organize files
#mv *.out $1
#mv *.dmp $1
mv *.tmp Out
mv *.aei Aei

### Print runtime
t2=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t2'))'

cd ..
dt=$(python -c $t2'-'$t1)
echo $1    $dt
