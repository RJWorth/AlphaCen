
### Location
dir='Prx02/Disk'	# or uncomment this to specify a directory
aeidir=paste('../Sims/Proxlike/',dir,'/Out/AeiOutFiles/',sep='')

### Generate list of disk particle names (for n particles)
n=100
disknames=rep( 'M', n)
	for (i in 1:n) disknames[i]=paste('M',i, sep='')
print(disknames)

### Determine dimensions
maxTlength = 0
for (obj in disknames)  {
	maxTlength=max(maxTlength,
					length(readLines(paste(aeidir,obj,'.aei',sep='') ) )-4) 
	}

ndims = dim(read.table(
				paste(aeidir,disknames[1],'.aei',sep=''), header=F,skip=4))[2]
if (ndims == 12) {
	cols=c('Time', 'a','e','i', 'mass','dens', 
									'x','y','z', 'vx','vy','vz')
	c.ind=c(2:4,7:12)
	} else {
	print('Wrong number of columns in .aei file')}

### Generate matrix of disk aei data
disk = array(data = NA, dim = c( n, maxTlength, ndims-3), 
	dimnames = list(NULL,NULL,cols[c.ind]) )
for (i in 1:length(disknames)) {
	Tlength = length(readLines( paste(aeidir,disknames[i],'.aei',sep='') ) )-4
	objdata = read.table(
		paste(aeidir,disknames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,c.ind]
	# It's stupid that I have to iterate over this, but otherwise it
	# turns my matrix into a list of n*maxTlength*ndim single elements...???
	for (j in 1:(ndims-3))	{
		for (k in 1:Tlength)	disk[i, k, j] = objdata[ k,j]	}
	# Get list of timesteps from the longest sim(s)
	if (Tlength>=maxTlength)	time=read.table(
		paste(aeidir,disknames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	}

sum3d = function(matrix3d)	{
	flattened = matrix3d[1,,]

	dims = dim(matrix3d)
	for (i in 2:dims[1])	{
		flattened = rbind(flattened,matrix3d[i,,])
	}
	
	summ = summary(flattened)
	return(summ)
	}

print(sum3d(disk))

extent=max(abs(disk[,,4:5]),na.rm=T)
grays = gray.colors(n, start = 0., end = 0.8)


pdf('DiskOrbits.pdf', height=10,width=10)
plot(0,0, pch=19, col='orange', 
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in 1:n) lines(disk[i,,4:5], col=grays[i])
#t=as.numeric(as.vector(stars[[1]]$Time))
dev.off()

