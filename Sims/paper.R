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
### Set up plot
marlist = list( c(0,1,1,0), c(0,0,1,1), c(1,1,0,0), c(1,0,0,1) )
xaxis = c(10, 1e5)
yaxis = c( 0, 1.0)

xaxlabs = c(10,100,1000,10000,100000)
yaxlabs = c(0,.2,.4,.6,.8,1)
### Side-by-side input and output parameters
##setEPS()
##postscript('../Presentation/ae_io.eps',width=5, height=5)
pdf('../Paper/Inserts/ae_io.pdf',width=6, height=6)
par(mfrow=c(2,2),oma=c( 2.75, 2.75, 2.5, 1.5))
options(scipen=5)

### Make 'before' plot
# make plot borders
if (version==1) par(mar=marlist[[1]]) else par(mar=marlist[[3]])

plot(.1,.1, pch='.', 
	log='x', axes=FALSE,frame.plot=TRUE,
	xlim=xaxis, ylim=yaxis, 
	xaxt='n',yaxt='n',
	main='',xlab='',ylab='')
if (version==1) {
	axis(1, tcl=0.25,cex.axis=.9,labels=FALSE)
	axis(2, tcl=0.25,at=yaxlabs,labels=as.character(yaxlabs))
	} else {	
	axis(1, tcl=0.25,cex.axis=.9,at=xaxlabs,labels=c(as.character(xaxlabs[-5]),''))
	axis(2, tcl=0.25,at=yaxlabs,labels=c(as.character(yaxlabs[-6]),''))
	}

# plot B where indB=='red'
#points(aB[indB=='red'], eB[indB=='red'], pch='.', col=indB[indB=='red'])
# plot B where indB isn't 'red'
points(aB[indB!='red'], eB[indB!='red'], pch='.', col=indB[indB!='red'])
# plot B from surviving systems (with black borders)
	for (j in 1:sum(surv))	{
	points(aB[surv][j], eB[surv][j], pch=20, col=allcols[1],cex=.6)
	points(aB[surv][j], eB[surv][j], pch=21, col='black',cex=.6)
	}
# plot all initial C parameters
points(aC, eC, pch='.', col=indC)
# plot C from surviving systems bigger
points(aC[surv], eC[surv], pch=20, col=indC[surv])
points(aC[prox], eC[prox], pch=21, col='orange')
# add something to proxima-like systems?
# lines connecting B and C from each surviving sim
	for (j in 1:sum(surv))	{
#	points(aB[prox][j], eB[prox][j], pch=21, col=colors()[257],cex=.7)
#	points(aC[prox][j], eC[prox][j], pch=21, col=colors()[257],cex=.9)
	lines(c(aB[surv][j], aC[surv][j]), c(eB[surv][j], eC[surv][j]),
	lty=3, col='black', lwd=.5)
	}

# add legend
if (version==1)	{
legend('bottomright', cex=0.9, ncol=1, 
	pch=c(20,20,15),
	col=c('red','blue','gray'),
	legend=c('Inner binary','Outer star','Proxima-like'))
	}
###############################################################################
### Make 'after' plot
if (version==1) par(mar=marlist[[2]]) else par(mar=marlist[[4]])

plot(.1,.1, 
	log='x', axes=FALSE,frame.plot=TRUE,
	xlim=xaxis, ylim=yaxis, 
	main='',xlab='',ylab='')

	axis(2, tcl=0.25,labels=FALSE)
if (version==2)	{
	axis(1, tcl=0.25,cex.axis=.9,at=xaxlabs,labels=as.character(xaxlabs))
	} else {
	axis(1, tcl=0.25,labels=FALSE)
	}

# Plot line for orbit with apocenter = 10,000 AU
polygon(c(a1.prx, a2.prx[length(a2.prx):1]), 
		c(e1.prx, e2.prx[length(a2.prx):1]), 
			col='gray',border='gray')
#lines( a1.prx, e1.prx, lty=2)
#lines( a2.prx[length(a2.prx):1], e2.prx[length(a2.prx):1], lty=2)
for (j in 1:sum(surv))	{
# Draw line from initial to final position
lines(c(aB[surv][j], aBf[surv][j]), c(eB[surv][j], eBf[surv][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[surv][j], aCf[surv][j]), c(eC[surv][j], eCf[surv][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(aBf[surv][j], aCf[surv][j]), c(eBf[surv][j], eCf[surv][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(surv))	{
	points(aBf[surv][j], eBf[surv][j], pch=20, col=allcols[1],cex=.6)
	points(aBf[surv][j], eBf[surv][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(aCf[surv], eCf[surv], pch=20, col=indC[surv])
points(aCf[prox], eCf[prox], pch=21, col='orange')

detach(SumAll)
	}

mtext(outer=TRUE,side=1,line= 1.5,cex=1.2,'Semimajor axis (AU)')
mtext(outer=TRUE,side=2,line= 1.3,cex=1.2,'Eccentricity')
mtext(outer=TRUE,side=3,line=0.75,cex=1.3,'Stellar Orbital Parameters')
mtext(outer=TRUE,side=3,line=-.75,cex=1.1,
	'Initial                              Final')
mtext(outer=TRUE,side=4,line=-.20,cex=1.1,
	'Equal masses               Current masses')

dev.off()

