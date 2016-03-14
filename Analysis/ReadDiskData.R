
print('Reading data')
### Generate list of disk particle names (for n particles)
print(dir)
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
	print('unrecognized name format => assuming disk is centered on A') }
### Determine dimensions
# find number of columns/variables
ndims = dim(read.table(
				paste(aeidir,disknames[1],'.aei',sep=''), header=F,skip=4))[2]
if (ndims == 12) {
	cols=c('Time', 'a','e','i', 'mass','dens', 
									'x','y','z', 'vx','vy','vz')
	} else {
	print('Wrong number of columns in .aei file')}
c.ind=c("a", "e", "i", "x", "y", "z", "vx","vy","vz")

# find longest time in this simulation
maxTlength = 0
time = array()
for (obj in c(disknames,starnames[-1]))  {
	thistime=read.table(
		paste(aeidir,obj,'.aei',sep=''), header=F,skip=4,
		col.names=cols)[,1]
	if (length(thistime)>length(time)) {
		time=thistime
		maxTlength=length(thistime)	} 
	}	# obj, disk/starnames

### Generate matrix of star aei data
star = array(data = 0., dim = c( nstars, maxTlength, length(c.ind)), 
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
	}
### Generate matrix of disk aei data
disk = ReadDisk(aeidir, disknames, n, maxTlength, c.ind)

print(sum3d(disk))


