### Get run version from input variable
args <- commandArgs(trailingOnly = F)
	print(length(args))
	print(args)
#if (length(args) >2) options(echo = FALSE)
readdir <- sub("-","",args[length(args)])

### If version not given at input, try using version specified here
if (length(dir) == 0 | 
	readdir=="-no-readline" | 
	readdir=="/usr/lib/R/bin/exec/R" | 
	readdir=="/usr/lib64/R/bin/exec/R" |
	readdir=="/usr/global/R/3.0.1/lib64/R/bin/exec/R" | 
	readdir=="/Library/Frameworks/R.framework/Resources/bin/exec/x86_64/R") {
		print('no args')
		readdir='Proxlike/Prx01/DiskB-2'
		}
	print(readdir)

### Locations
afile = paste(readdir,'/Out/AeiOutFiles/a.out',sep='')
efile = paste(readdir,'/Out/AeiOutFiles/e.out',sep='')
#ifile = paste(readdir,'/Out/AeiOutFiles/i.out',sep='')

xfile = paste(readdir,'/Out/AeiOutFiles/x.out',sep='')
yfile = paste(readdir,'/Out/AeiOutFiles/y.out',sep='')
zfile = paste(readdir,'/Out/AeiOutFiles/z.out',sep='')

tfile = paste(readdir,'/Out/AeiOutFiles/t.out',sep='')

### Read in data
t=read.table(tfile,header=F)[[1]]
a=read.table(afile,header=T)
e=read.table(efile,header=T)
x=read.table(xfile,header=T)
y=read.table(yfile,header=T)
z=read.table(zfile,header=T)

### Get object names from column names
objs=colnames(a)
n = length(grep('M',objs))

### Calculate other quantities
nt=dim(a)[1]						# number of timeesteps

last = max(which(a[nt, -1:-3]>0))+3	# index of outermost surviving disk obj
aout = a[nt,last]					# a of last obj
etop = max(e[,4:last])				# max e of objs inside last ind

p=a*(1-e)							# perihelion
r=sqrt(x^2+y^2+z^2)					# dist (from A?)

### Calculate gas disk mass over time (in mSun)
Mg = 10*(1.-t/(440000))
for (i in 1:length(Mg)) Mg[i] = max(Mg[i],0.)


#plot(a$PrxCen,e$PrxCen,type='l')
if (grep('-3',readdir)==1) {
	source('../Analysis/DiskAEPlot3.R')
} else {
	source('../Analysis/DiskAEPlot2.R')

