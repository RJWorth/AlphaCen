
maxtime=10
mintime=5
seq=$(jot $(echo "$maxtime-$mintime+1" | bc) $mintime)
echo $seq
