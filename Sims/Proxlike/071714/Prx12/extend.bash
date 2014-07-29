#!/bin/bash
############################################################################### 
### after modifying param.dmp, run this to restart the simulation
### with folder name as first argument, merc name as second. Example:
### ./extend.bash 4-1 CDir10 > run.pipe &

t1=$(date +%s)

cd $1/Out
	./merc_AC$2
	./elem
	\mv *.aei AeiOutFiles
cd ../..

t2=$(date +%s)
echo $1"	"$(echo "$t2 - $t1"|bc ) >> runtime.txt

cd ../..; ./email.sh $1 '' 'AC extension finished'; cd Proxlike

