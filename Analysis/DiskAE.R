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
		readdir='Proxlike/Prx01/DiskB-3'
		}
	print(readdir)

### Locations
afile = paste(readdir,'/Out/AeiOutFiles/a.out',sep='')
efile = paste(readdir,'/Out/AeiOutFiles/e.out',sep='')
#ifile = paste(readdir,'/Out/AeiOutFiles/i.out',sep='')

xfile = paste(readdir,'/Out/AeiOutFiles/x.out',sep='')
yfile = paste(readdir,'/Out/AeiOutFiles/y.out',sep='')
zfile = paste(readdir,'/Out/AeiOutFiles/z.out',sep='')


a=read.table(afile,header=T)
e=read.table(efile,header=T)
x=read.table(xfile,header=T)
y=read.table(yfile,header=T)
z=read.table(zfile,header=T)
nt=dim(a)[1]

p=a*(1-e)
r=sqrt(x^2+y^2+z^2)

#plot(a$PrxCen,e$PrxCen,type='l')
source('../Analysis/DiskAEPlot.R')
