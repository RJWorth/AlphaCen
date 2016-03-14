### Include file to use with M/R-BSDisk.bash scripts

### Run where?
dir1='BS/Prx'
dir2='/Disk'
DirNums=(17 21 28 34)
echo $dir1 ${DirNums[*]}

### How clean?
newrun=T

newdisk=F   # generate a new disk, T or F
#amin=0.1   # minimum extent of disk, in AU
#sz=small       # add to big.in or small.in? (not used yet)

newmerc=T   # recompile the merc/elem executables?
vers='ury_TG.for'   # 'merc'+vers = filename for mercury

### Simulation parameters
mA=0.123    #1.105
mB=0.123    #0.934
mC=0.123

mintime=3   # = log(years)
maxtime=7   # = log(years)
output=1    # = log(years)
step=0.0003 # = days
user='yes'  # use user-defined forces?
alg='bs'    # integrator algorithm

### Range for iterations
machine=$(hostname -s)
if [ $machine = chloe ]; then
    timerange=$(jot $(echo "$maxtime-$mintime+1" | bc) $mintime)
else
    timerange=$(seq $mintime $maxtime)
fi
echo '  timerange '$timerange

