#!/bin/bash
set -e  # stop script on error
############################################################################### 

### Sim parameters
incl=master.inc
source $incl

### Determine relevant directories
if [ $machine = Mirka ]; then
Dirs=()
elif [ $machine = chloe ];then
Dirs=(C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12)
Dirs=()
elif [ $machine = shapiro ]; then
Dirs=(S01 S02 S03 S04 S05 S06 S07)
Dirs=()
else
echo '*** unrecognized machine; not running anything ***'
fi 

### Start run for each directory on this machine
echo $Dirs

for i in ${Dirs[*]}; do	
	nice -n 10 ./run.bash $i $incl > $i/run.pipe 2>&1 &
	echo 'started '$i': '$!
done	# i in Dirs

