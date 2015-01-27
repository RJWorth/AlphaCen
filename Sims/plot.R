### Make plots and such for the suite of simulations indicated.
### run from command line:
#   R CMD BATCH -1 plot.R

### Get run version from input variable
args <- commandArgs(trailingOnly = F)
	print(length(args))
	print(args)
if (length(args) >2) options(echo = FALSE)
version <- sub("-","",args[length(args)])
	print(version)

### If version not specified, try using version specified here
if (length(version) == 0 | 
	version=="-no-readline" | 
	version=="/usr/lib/R/bin/exec/R" | 
	version=="/usr/lib64/R/bin/exec/R") {
		print('no args')
		version = 1	}

### Packages
library(xtable)

### Settings
r.prx 	= 10000.		# min apoapse to be called 'proxima-like'
r.big 	= 20000.		# max apoapse to be called 'proxima-like'
DoCalcRunTime = FALSE	# Whether to run the calc time script

### Files to read data from
### Regular (current) version:
if (version==1)	{
	nmissing = 10				# number of sims lost on the 1e9 step
	case  = 'Equal masses'
	prefix='../Paper/Inserts/EqualMasses/'	# Prefix for plot files
	sumfiles = c('SumAll.out')
	} else if (version==2)	{
### Older versions of sims (eps instead of E, no dE, unequal masses)
	nmissing = 0				# number of sims lost on the 1e9 step
	case  = 'True masses'
	prefix='../Paper/Inserts/TrueMasses/'	# Prefix for plot files
	sumfiles=c(	'Saved/CurrentMasses/SumAll081414.out',
				'Saved/CurrentMasses/SumAll080614.out',
				'Saved/CurrentMasses/SumAll072314.out' )
	} else print('invalid version')

###############################################################################
### Calculate average runtime per sim
if (DoCalcRunTime) source('CalcRunTime.R')

###############################################################################
### Read in Summary.out, normal
source('ReadSumFiles.R')
nsims=dim(SumAll)[1]

###############################################################################
### Calculate various plotting variables and such
source('PlotVars.R')

###############################################################################
### Create all the pdf plots
source('MakePlots.R')

###############################################################################
### Print proxima-like systems to a file (aCf>acutoff)
#acutoff=5000	# AU
#sink(paste(prefix,'proxlike.txt',sep=''))
#options(width=300)
#print(SumAll[prox,],row.names=FALSE)
#options(width=80)
#sink()

###############################################################################
###################### Print tables for latex file
### Print summary of these simulations by fate of system
### read in current written version
SumTableFile='../Paper/Inserts/SumTable.tex'
SumTable=read.table(SumTableFile, colClasses='character',
	sep='&',strip.white=TRUE, 
	skip=2, comment.char='\\',
	row.names=1,header=TRUE, check.names=FALSE)
### update data in row for this version #
#SumTable$'Broken'       =rep( m$brkn, 2)
SumTable[,case]=c(toString(nsims+nmissing), 
	toString( m$surv),toString( m$grow),toString(m$prox),toString(m$huge), 	
	toString(m$singB),toString(m$singC),toString(m$doub),toString(m$coll),
	toString(m$brkn))

xSumTable=xtable(SumTable)
### write updated table back to the file
print(xSumTable, file=SumTableFile,
	only.contents    = TRUE,
#	include.colnames = FALSE,
	hline.after      = c(0,2),
	  )
write('',file=SumTableFile,append=TRUE)	# adds EOF
#------------------------------------------------------------------------------
### Print latex table of binary parameters in proxima-like systems

Bparams = data.frame( rep(' ',length(aB)),aB,eB,iB, aBf, eBf, iBf)[prox1,]
	colnames(Bparams)[1] = ' '
print(xtable(Bparams[order(Bparams$aB),]), 
	file=paste(prefix,'Bparams.tex',sep=''),
	only.contents    = TRUE,
	include.rownames = FALSE,
	include.colnames = FALSE,
	hline.after      = NULL  )

Cparams = data.frame( rep(' ',length(aB)), aC,eC,iC, aCf, eCf, iCf)[prox1,]
	colnames(Cparams)[1] = ' '
print(xtable(Cparams[order(Bparams$aB),]), 
	file=paste(prefix,'Cparams.tex',sep=''),
	only.contents    = TRUE,
	include.rownames = FALSE,
	include.colnames = FALSE,
	hline.after      = NULL  )

###############################################################################
detach(SumAll)

### If called from shell command line, force exit from R
if (length(args) > 2) q('no')


