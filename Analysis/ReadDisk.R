
### Location
dir='Prx02/Disk'	# or uncomment this to specify a directory
aeidir=paste('../Sims/Proxlike/',dir,'/Out/AeiOutFiles/',sep='')

### Load relevant constants/functions etc
source('../Analysis/DiskUtils.R')

### Generate list of disk particle names (for n particles)
n=100
disknames=rep( 'M', n)
	for (i in 1:n) disknames[i]=paste('M',i, sep='')
print(disknames)

starnames=c('AlCenA','AlCenB')
nstars=length(starnames)

### Determine dimensions
maxTlength = 0
for (obj in c(disknames,starnames[-1]))  {
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

### Generate matrix of star aei data
star = array(data = 0., dim = c( nstars, maxTlength, ndims-3), 
	dimnames = list(NULL,NULL,cols[c.ind]) )
for (i in 2:length(starnames)) {
	Tlength = length(readLines( paste(aeidir,starnames[i],'.aei',sep='') ) )-4
	objdata = read.table(
		paste(aeidir,starnames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,c.ind]
	# It's stupid that I have to iterate over this, but otherwise it
	# turns my matrix into a list of n*maxTlength*ndim single elements...???
	for (j in 1:(ndims-3))	{
		for (k in 1:Tlength)	star[i, k, j] = objdata[ k,j]	}
	# Get list of timesteps from the longest sim(s)
	if (Tlength>length(time))	time=read.table(
		paste(aeidir,starnames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	}
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
	if (Tlength>length(time))	time=read.table(
		paste(aeidir,disknames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	}

print(sum3d(disk))

#extent=max(abs(c(disk[,,4:5],star[,,4:5])),na.rm=T)
extent=max(abs(c(disk[,1,4:5],star[,,4:5])),na.rm=T)
grays = gray.colors(n, start = 0., end = 0.8)

### Measure the difference between actual and expected velocity, 
### as a fraction of expected
verr=(v(1:100,1)-vorb_expect(1:100,1))/vorb_expect(1:100,1)

### Make a 2D array of initial parameters for each object
pairframe=cbind(disk[,1,],verr)

### Make plots
pdf('DiskOrbits_Both.pdf', height=10,width=10)
pairs(pairframe,pch=20,col=grays)

#plot(r(1,1:1001), keps(1,1:1001), pch=20)

plot(0,0, pch=19, col='orange', 
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in 1:nstars) lines(star[i,,4:5], col='orange')
for (i in 1:n)      lines(disk[i,,4:5], col=grays[i])
#t=as.numeric(as.vector(stars[[1]]$Time))
dev.off()

