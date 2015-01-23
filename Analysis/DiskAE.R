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
		readdir='Proxlike/Prx01'	#'Prx02/Disk'
		}
	print(readdir)

### Locations
if (substr(readdir,1,1)=='/') simdir=readdir else {
	simdir=paste('../Sims/',readdir,sep='') }

#subdirs=c('/DiskA-2','/DiskA-3','/DiskB-2','/DiskB-3')
subdirs=c('/DiskB-3')	#'/DiskB-2','/DiskB-3')

### Load relevant constants/functions etc
source('../Analysis/DiskUtils.R')
source('../Analysis/ReadDiskFns.R')

#survgrid = array(data = NA, dim = c(length(subdirs),4), 
#	dimnames = list(subdirs,c('surv','stab','rEdge','rAbin')) )

for (dirnum in 1:length(subdirs))	{
	dir=paste(readdir,subdirs[dirnum],sep='')
	aeidir=paste(simdir,subdirs[dirnum],'/Out/AeiOutFiles/',sep='')
	print(aeidir)

	source('../Analysis/ReadDiskData.R')
	}





