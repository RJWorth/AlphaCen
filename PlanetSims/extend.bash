#!/bin/bash
set -e  # stop script on error
##############################################################################
# This script extends the most recent Merc95 by a single stage
# ./extend.bash dir tf >> dir/run.pipe &

t1=$(date +%s)

machine=$(hostname -s)
dir=$1

tf=$2	# = log(years)

mStar=0.934
rEj=1e2

home=..
c=../Code
out=Out/

### Move to working dir
cd $dir

echo '==================== extend t = 1e'$tf' yrs ===================='

### Write param.dmp file
tOut=$(python -c 'print('$tf'-3)')
echo $tOut
exit 1
python -c "import Merc as M; M.WriteParamInFile(loc='"$out"', f='dmp', tf=365.25e"$tf", tOut=365.25e"$tOut", mStar="$mStar", rEj="$rEj")"

#### Run mercury
./merc_$1
echo '----------------------- t = 1e'$tf' yrs -------------------------'

### Print runtime thus far
t2=$(date +%s)
dt=$(python -c 'print('$t2'-'$t1')')
echo $dir'    '$tf'    '$dt >> $home/runtimes.txt

echo '====================== end t = 1e'$tf' yrs ===================='


### Run element
./elem_$1

### Organize files
mv *.tmp Out
mv *.aei Aei

### Print runtime
t3=$(date +%s)
python -c 'import Merc; print(Merc.WriteRuntime('$t1','$t3'))'

cd $home

dt=$(python -c 'print('$t3'-'$t1')')
echo $dir'    '$tf'    '$dt >> runtimes.txt


