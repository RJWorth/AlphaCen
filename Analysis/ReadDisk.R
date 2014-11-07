### Get run version from input variable
args <- commandArgs(trailingOnly = F)
	print(length(args))
	print(args)
if (length(args) >2) options(echo = FALSE)
dir <- sub("-","",args[length(args)])

### If version not given at input, try using version specified here
if (length(dir) == 0 | 
	dir=="-no-readline" | 
	dir=="/usr/lib/R/bin/exec/R" | 
	dir=="/usr/lib64/R/bin/exec/R" |
	dir=="/Library/Frameworks/R.framework/Resources/bin/exec/x86_64/R") {
		print('no args')
		dir='Proxlike/Prx01/DiskB-2'	#'Prx02/Disk'
		}
	print(dir)

### Locations
if (substr(dir,1,1)=='/') simdir=dir else {
	simdir=paste('../Sims/',dir,sep='') }
aeidir=paste(simdir,'/Out/AeiOutFiles/',sep='')

### Load relevant constants/functions etc
source('../Analysis/DiskUtils.R')

### Generate list of disk particle names (for n particles)
n=(length(readLines(paste(dir,'/In/small.in',sep='')))-5)/4
disknames=rep( 'M', n)
	if (dir == 'Proxlike/071714/Prx1/Disk4')	{
		for (i in 1:n) disknames[i]=paste('M',i+100, sep='') } else {
		for (i in 1:n) disknames[i]=paste('M',i, sep='') }
print(disknames)

starnames=c('AlCenA','AlCenB')
nstars=length(starnames)

### Determine where disk is centered
cent = substr(dir, nchar(dir)-2, nchar(dir)-2)
if (cent == 'A') cent=1 else if (cent =='B') cent=2 else {
	print('unrecognized name format, assuming disk is centered on A') }
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
	} else {
	print('Wrong number of columns in .aei file')}
c.ind=c("a", "e", "i", "x", "y", "z", "vx","vy","vz")

### Generate matrix of star aei data
star = array(data = NA, dim = c( nstars, maxTlength, length(c.ind)), 
	dimnames = list(NULL,NULL,c.ind) )
for (i in 2:length(starnames)) {
	Tlength = length(readLines( paste(aeidir,starnames[i],'.aei',sep='') ) )-4
	objdata = read.table(
		paste(aeidir,starnames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,c.ind]
	# It's stupid that I have to iterate over this, but otherwise it
	# turns my matrix into a list of n*maxTlength*ndim single elements...???
	for (j in 1:length(c.ind))	{
		for (k in 1:Tlength)	star[i, k, j] = objdata[ k,j]	}
	# Get list of timesteps from the longest sim(s)
	if (Tlength>length(time))	time=read.table(
		paste(aeidir,starnames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	}
### Generate matrix of disk aei data
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
			}
		r = sqrt( (disk[i,k,'x']-star[cent,k,'x'])^2 +
				  (disk[i,k,'y']-star[cent,k,'y'])^2 +
				  (disk[i,k,'z']-star[cent,k,'z'])^2 )
		disk[i, k, length(c.ind)+1] = r
	}
	# Get list of timesteps from the longest sim(s)
	if (Tlength>length(time))	time=read.table(
		paste(aeidir,disknames[i],'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	}

print(sum3d(disk))

### Calculate % of disk surviving, as a f'n of time
outerbound = min(star[2,,'a'],na.rm=T)
surviving = matrix(, nrow = length(time), ncol = dim(disk)[1])
stable = matrix(, nrow = length(time), ncol = dim(disk)[1])
surv.per = rep(0., length(time))
stab.per = rep(0., length(time))
for (i in 1:length(disknames)) {
	surviving[,i] =  !is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
	stable[,i]    = (!is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
					 & (abs(disk[i,,'r']-disk[i,1,'r']) <= 1.)
					 & ((cent==2) | (disk[i,,'a'] >= 0.)))
	}
for (j in 1:length(time))	{
	surv.per[j] = sum(surviving[j,])/n
	stab.per[j] = sum(   stable[j,])/n
	}
#extent=max(abs(c(disk[,,4:5],star[,,4:5])),na.rm=T)
extent=max(abs(c(disk[,1,4:5],star[,,4:5])),na.rm=T)
grays = gray.colors(n, start = 0., end = 0.8)
heat  = heat.colors(n)

### Resample log scale for time
samprate=100
logt = min(log10(time[-1])) + 
		(max(log10(time))-min(log10(time[-1])))*(0:(samprate-1))/(samprate-1)
resampt = 10^logt
### OR linear scale
#resampt = min(time) + 
#			(max(time)-min(time))*(0:(samprate-1))/(samprate-1)


### Make image matrix of disk survival
diskimg = matrix(, nrow = samprate, ncol = dim(disk)[1])
	for (j in 1:n)	{
		for (k in 1:samprate)	{
			i = which( abs(time-resampt[k]) == min(abs(time-resampt[k])) )
			if (length(i) > 1) i=min(i)
			if ((surviving[i,j]==TRUE) & (stable[i,j] == TRUE)) {
				diskimg[k,j]=j
			} else if ((surviving[i,j]==TRUE) & (stable[i,j] == FALSE)) {
				diskimg[k,j]=1+n
			} else if ((surviving[i,j]==FALSE)) {
				diskimg[k,j]=2+n
			} else {print(paste(k,j))}
	}}

### Measure the difference between actual and expected velocity, 
### as a fraction of expected
#verr=(v(1:100,1)-vorb_expect(1:100,1))/vorb_expect(1:100,1)

### Make a 2D array of initial parameters for each object
#pairframe=cbind(disk[,1,],verr)
pairframe=disk[,1,]

time.all  =               1:maxTlength
timeslice = (maxTlength-99):maxTlength
### Make plots
pdf(paste(simdir,'/DiskPairs.pdf',sep=''), height=10,width=10)
pairs(pairframe,pch=20,col=grays)
dev.off()
###
pdf(paste(simdir,'/DiskSurv.pdf',sep=''), height=4,width=8)

image(diskimg,col=c(heat,'lightgray','darkgray'), axes=FALSE)
axis(4, col="red", lwd=2,
	labels=c( min(disk[,1,'r']), 5, 10 ),
	at=c(0,.5,1) )

par(new=T)
plot(time[-1],surv.per[-1], type='l',lwd=2,col='blue',log='x',
	xaxs='i',yaxs='i',
	xlim=c(min(time[-1]),max(time)),ylim=c(0,1),
	ylab='Surviving percentage',xlab='Time (yrs)',
	main=paste(stab.per[length(time)]*100,'% of disk remains stable',sep=''))
#lines(time,surv.per, lwd=2)
lines(time,stab.per, lwd=3, lty=3,col='blue')
legend('bottomleft',
	legend=c('Fraction surviving','Fraction remaining stable'),
	lty=c(1,3),col='blue',lwd=3)

dev.off()
###
pdf(paste(simdir,'/DiskEvol.pdf',sep=''), height=15,width=15)
par(mfrow=c(3,1))

for (column in c('x','y','z'))	{
plot(time[-1],star[2,time.all[-1],column], type='l', col='orange',log='x',
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(time[-1],disk[i,time.all[-1],column], col=grays[i])
	}

dev.off()
###
pdf(paste(simdir,'/DiskOrbits.pdf',sep=''), height=10,width=10)

plot(0,0, pch=19, col='orange', main='End of sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      points(disk[i,timeslice,4:5], col=grays[i],pch=20)
for (i in 2:nstars) points(star[i,timeslice,4:5], col='orange',pch=20)

plot(0,0, pch=19, col='orange', main='Whole sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(disk[i, time.all,4:5], col=grays[i])
for (i in 2:nstars) lines(star[i, time.all,4:5], col='orange')

plot(0,0, pch=19, col='orange', main='Start of sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(disk[i,     1:50,4:5], col=grays[i])
for (i in 2:nstars) lines(star[i,     1:50,4:5], col='orange')
dev.off()

### If called from shell command line, force exit from R
if (length(args) > 2) q('no')

