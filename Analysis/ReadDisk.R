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
		print('No args given: using defaut directory!')
		readdir='Proxlike/072314/Prx24'
		}
	print(readdir)

### Locations
if (substr(readdir,1,1)=='/') simdir=readdir else {
	simdir=paste('../Sims/',readdir,sep='') }

#subdirs=c('/DiskA-2','/DiskA-3','/DiskB-2','/DiskB-3')
subdirs=c('/DiskB-2','/DiskB-3')
#subdirs=c('/DiskA-2','/DiskA-3')

### Load relevant constants/functions etc
source('../Analysis/DiskUtils.R')
source('../Analysis/ReadDiskFns.R')

survgrid = array(data = NA, dim = c(length(subdirs),4), 
	dimnames = list(subdirs,c('surv','stab','rEdge','rAbin')) )

l.disk     = list()
l.diskimg  = list()
l.time     = list()
l.surv.per = list()
l.stab.per = list()
l.edge     = list()

for (dirnum in 1:length(subdirs))	{
	dir=paste(readdir,subdirs[dirnum],sep='')
	aeidir=paste(simdir,subdirs[dirnum],'/Out/AeiOutFiles/',sep='')
	print(aeidir)

	source('../Analysis/ReadDiskData.R')
	source('../Analysis/ReadDiskAnalysis.R')
#	source('../Analysis/ReadDiskPlots.R')

	l.disk[[dirnum]]     = disk
	l.diskimg[[dirnum]]  = diskimg
	l.time[[dirnum]]     = time
	l.surv.per[[dirnum]] = surv.per
	l.stab.per[[dirnum]] = stab.per
	l.edge[[dirnum]]     = edge
	
	survgrid[dirnum,1] = l.surv.per[[dirnum]][length(l.stab.per[[dirnum]])]
	survgrid[dirnum,2] = l.stab.per[[dirnum]][length(l.stab.per[[dirnum]])]
	survgrid[dirnum,3] = r0[edge[length(edge)]]
	survgrid[dirnum,4] = rA[edge[length(edge)]]

	} # dirnum in subdirs
sink(paste(simdir,'/SurvGrid.txt',sep=''))
print(survgrid)
sink()

source('../Analysis/ReadDiskPlot4.R')

################################## Cleanup ####################################
### If called from shell command line, force exit from R
if (length(args) > 2) q('no')

