ReadDisk = function(aeidir, disknames, n, maxTlength, c.ind)	{
disk = array(data = NA, dim = c( n, maxTlength, length(c.ind)+1), 
	dimnames = list(NULL,NULL,c(c.ind,'r')) )
for (i in 1:length(disknames)) {
	Tlength = length(readLines( paste(aeidir,disknames[i],'.aei',sep='') ) )-4
	objdata = read.table(
		paste(aeidir,disknames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,c.ind]
	# It's stupid that I have to iterate over this, but otherwise it
	# turns my matrix into a list of n*maxTlength*ndim single elements...???
	for (k in 1:Tlength)	{
		for (j in 1:length(c.ind))	{
				disk[i, k, j] = objdata[k,j]
			}	# j, c.ind
		r = sqrt( (disk[i,k,'x']-star[cent,k,'x'])^2 +
				  (disk[i,k,'y']-star[cent,k,'y'])^2 +
				  (disk[i,k,'z']-star[cent,k,'z'])^2 )
		disk[i, k, length(c.ind)+1] = r
	}	# k, Tlength
	}	# i, disknames
	return(disk)
	}	# function
#------------------------------------------------------------------------------
FindClosest = function(x,n)	{
	Delta  = abs(x-n)
	MinDel = min(Delta)
	MinInd = which(Delta == MinDel)
	if (length(MinInd) >1) MinInd = MinInd[1]
	return(MinInd)
	}	# function
#------------------------------------------------------------------------------
### Moving average that gives NAs for boundaries
#MovingAvg <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
#------------------------------------------------------------------------------
### Moving average that gives one-sided average for boundaries.
### n = number of points on either side to include (so n=1 avgs 3 points)
### (needs work for n>1 -- will fail)
MovingAvg = function(x,n=1)	{
	len=length(x)
	avg = rep(NA,len)
	for (i in 1:len)	{
		if (i == 1) avg[i] = mean(x[i:(i+n)]) else if (i == len) {
			avg[i]=mean(x[(i-1):i])	} else {
			avg[i] = mean(x[(i-n):(i+n)])	}
		} 
	return(avg)
	}

