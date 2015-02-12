### Make plots and such for the suite of simulations indicated.
### run from command line:
#   R CMD BATCH poster.R

### Packages
library(xtable)

### Settings
r.prx 	= 10000.		# min apoapse to be called 'proxima-like'
r.big 	= 20000.		# max apoapse to be called 'proxima-like'
DoCalcRunTime = FALSE	# Whether to run the calc time script


###############################################################################
### Read in Summary.out for each case
nmissing=c(0,10)

### Late-interaction, true (high) mass case
#	nmissing = 0				# number of sims lost on the 1e9 step
#	case  = 'True masses'
	prefix='../Presentation/True'	# Prefix for plot files
	sumfiles=c(	'Saved/CurrentMasses/SumAll081414.out',
				'Saved/CurrentMasses/SumAll080614.out',
				'Saved/CurrentMasses/SumAll072314.out' )
source('ReadSumFiles.R')
nsims.l=dim(SumAll)[1]
Case=rep('Late',nsims.l)
detach(SumAll)
StarData = cbind(Case,SumAll)

### Early-interaction, equal (low) mass case
#	nmissing = 10				# number of sims lost on the 1e9 step
#	case  = 'Equal masses'
	prefix='../Presentation/Eql'	# Prefix for plot files
	sumfiles = c('SumAll.out') 
source('ReadSumFiles.R')
nsims.e=dim(SumAll)[1]
Case=rep('Early',nsims.e)
detach(SumAll)

nsims = c(nsims.l,nsims.e)

StarData = rbind(StarData, cbind(Case,SumAll[,c(-22,-23)]) )

pB  = StarData$aB*(1-StarData$eB)
apB = StarData$aB*(1+StarData$eB)
pC  = StarData$aC*(1-StarData$eC)
apC = StarData$aC*(1+StarData$eC)

pBf  = StarData$aBf*(1-StarData$eBf)
apBf = StarData$aBf*(1+StarData$eBf)
pCf  = StarData$aCf*(1-StarData$eCf)
apCf = StarData$aCf*(1+StarData$eCf)

StarData = cbind(StarData, pB,apB,pC,apC, pBf,apBf,pCf,apCf)

rm(Case)
attach(StarData)

###############################################################################
### Calculate various plotting variables and such
source('PlotVars.R')

###############################################################################
### Set up stellar plot
source('PaperIOPlot.R')

###############################################################################
### Read in disk data (for 03?)
#source('../Analysis/ReadDisk.R')

###############################################################################







