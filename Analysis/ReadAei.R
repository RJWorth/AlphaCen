#read .aei for stars
print('Read .aei files')

args <- commandArgs(trailingOnly = F)
dir  <- sub("-","",args[length(args)])

### If version not given at input, try using version specified here
if (length(dir) == 0 | 
	dir=="-no-readline" | 
	dir=="/usr/lib/R/bin/exec/R" | 
	dir=="/usr/lib64/R/bin/exec/R" |
	dir=="/Library/Frameworks/R.framework/Resources/bin/exec/x86_64/R") {
		print('no args')
		dir='Proxlike/Prx04/Original'
		}
print(dir)

#aeidir=paste('../Sims/Proxlike/',dir,'/Out/AeiOutFiles/',sep='')
aeidir=paste('../Sims/',dir,'/Out/AeiOutFiles/',sep='')	#location of directory

mode = 'triple'	# binary or triple
### Read aei for stars
if (mode=='triple') starnames=c('AlCenB','PrxCen') else
	if (mode=='binary') starnames=c('AlCenB')
stars=list()
for (j in 1:length(starnames)) stars[[j]]=read.table(
	paste(aeidir,starnames[j],'.aei',sep=''), header=F,skip=4,
	col.names=c('Time','a','e','i','mass','dens', 'x','y','z','vx','vy','vz')
	)[,c(1:3,7:12)]
tB=as.numeric(as.vector(stars[[1]]$Time))
tC=as.numeric(as.vector(stars[[2]]$Time))
if (length(tB) >= length(tC)) t = tB else t = tC

#disknames=rep( 'M', 100)
#	for (i in 1:100) disknames[i]=paste('M',i, sep='')
#print(disknames)
#disk=list()
#for (j in 1:length(disknames)) disk[[j]]=read.table(
#	paste(aeidir,disknames[j],'.aei',sep=''), header=F,skip=4,
#	col.names=c('Time','a','e','i','mass','dens', 'x','y','z','vx','vy','vz')
#	)[,c(1:3,7:12)]
#print(summary(disk))

### Get x and v in center-of-momentum frame
### Define masses
mSun	= 1.9891e30		# kg
if (mode=='triple')		m = c(1.105*mSun, 0.934*mSun, 0.123*mSun) else 
	if (mode=='binary')	m = c(1.105*mSun, 0.934*mSun, 0.)
M = sum(m)

### Find CM of AB and ABC
CM2=array(dim=c(length(t),6))
CM3=array(dim=c(length(t),6))
for (k in 1:length(t))	{
	for (l in (1:6))	{
		CM2[k,l]=(1/M)*(m[2]*stars[[1]][k,l+3])
		if (mode=='triple') {
			CM3[k,l]=(1/M)*(m[2]*stars[[1]][k,l+3] +
						    m[3]*stars[[2]][k,l+3] )}	}}

### Get x and v relative to CM
stars2=stars
stars3=stars
for (j in 1:length(starnames)) {	
	for (k in 1:length(t))	{
		stars2[[j]][k,4:9]=stars[[j]][k,4:9]-CM2[k,]
		if (mode=='triple') stars3[[j]][k,4:9]=stars[[j]][k,4:9]-CM3[k,]
		}}

#print(CM[1,])
#print(stars[[1]][7401,])
#print(stars[[1]]$x[7401])

### save number of timesteps in file to read later
sink(file='lengthT.txt')
cat(length(t))
sink()
### Reduce number of timesteps
redux=10
#DoAEPlot=(1:length(t))[ ((1:length(t))%%redux)==0]
# extension indices
#DoAEPlot=(6401:7401)[(1:1000%%1)==0]
start=5-4		#start=6405-4
stop =5505-4	#stop=7405-4
range=start:stop
DoAEPlot=range[(1:length(range)%%redux)==0]

sink(file='steps.txt')
cat(DoAEPlot[1])	# write first element to one line
cat('\n')
cat(DoAEPlot[-1])	# write the rest of the elements to the next line
sink()


### x-y plots
if ( length(DoAEPlot)>1 | DoAEPlot[1]>0 ) source('AEPlot.R')

### If called from shell command line, force exit from R
if (length(args) > 2) q('no')


