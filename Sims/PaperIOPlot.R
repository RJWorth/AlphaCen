
marlist = list( c(0,1,1,0), c(0,0,1,1), c(1,1,0,0), c(1,0,0,1) )
xaxis = c(10, 1e5)
yaxis = c( 0, 1.0)

xaxlabs = c(10,100,1000,10000,100000)
yaxlabs = c(0,.2,.4,.6,.8,1)

### Side-by-side input and output parameters
##setEPS()
##postscript('../Paper/Inserts/ae_io.eps',width=6, height=6)
pdf('../Paper/Inserts/ae_io.pdf',width=6, height=6)
par(mfrow=c(2,2),oma=c( 2.75, 2.75, 2.5, 1.5))
options(scipen=5)

### Make 'before' plot
# make plot borders
for (k in 1:2)	{
c= Case==levels(Case)[k]

if (k==1) par(mar=marlist[[1]]) else par(mar=marlist[[3]])

plot(.1,.1, pch='.', 
	log='x', axes=FALSE,frame.plot=TRUE,
	xlim=xaxis, ylim=yaxis, 
	xaxt='n',yaxt='n',
	main='',xlab='',ylab='')
if (k==1) {
	axis(1, tcl=0.25,cex.axis=.9,labels=FALSE)
	axis(2, tcl=0.25,at=yaxlabs,labels=as.character(yaxlabs))
	} else {	
	axis(1, tcl=0.25,cex.axis=.9,at=xaxlabs,labels=c(as.character(xaxlabs[-5]),''))
	axis(2, tcl=0.25,at=yaxlabs,labels=c(as.character(yaxlabs[-6]),''))
	}

# plot B where indB=='red'
#points(aB[indB=='red'], eB[indB=='red'], pch='.', col=indB[indB=='red'])
# plot B where indB isn't 'red'
points(aB[indB!='red' & c], 
	   eB[indB!='red' & c], pch='.', col=allcols[1])
# plot B from surviving systems (with black borders)
	for (j in 1:sum(surv[c]))	{
	points(aB[surv & c][j], eB[surv & c][j], pch=20, col=allcols[1],cex=.6)
	points(aB[surv & c][j], eB[surv & c][j], pch=21, col='black',cex=.6)
	}	# end j in 1:sum(surv[c])
# plot all initial C parameters
points(aC[c], eC[c], pch='.', col=allcols[4])
# plot C from surviving systems bigger
points(aC[surv & c], eC[surv & c], pch=20, col=allcols[4])
points(aC[prox & c], eC[prox & c], pch=21, col=allcols[2])
# add something to proxima-like systems?
# lines connecting B and C from each surviving sim
	for (j in 1:sum(surv[c]))	{
#	points(aB[prox][j], eB[prox][j], pch=21, col=colors()[257],cex=.7)
#	points(aC[prox][j], eC[prox][j], pch=21, col=colors()[257],cex=.9)
	lines(c(aB[surv & c][j], aC[surv & c][j]), 
		  c(eB[surv & c][j], eC[surv & c][j]),
		  lty=3, col='black', lwd=.5)
	}	# end j in 1:sum(surv[c])

# add legend
if (k==1)	{
legend('bottomright', cex=0.9, ncol=1, 
	pch=c(20,20,15),
	col=c(allcols[1],allcols[4],'gray'),
	legend=c('Inner binary','Outer star','Proxima-like'))
	}
###############################################################################
### Make 'after' plot
if (k==1) par(mar=marlist[[2]]) else par(mar=marlist[[4]])

plot(.1,.1, 
	log='x', axes=FALSE,frame.plot=TRUE,
	xlim=xaxis, ylim=yaxis, 
	main='',xlab='',ylab='')

	axis(2, tcl=0.25,labels=FALSE)
if (k==2)	{
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
for (j in 1:sum(surv[c]))	{
# Draw line from initial to final position
lines(c(aB[surv & c][j], aBf[surv & c][j]), c(eB[surv & c][j], eBf[surv & c][j]),
	col=allcols[2], lwd=.5)
lines(c(aC[surv & c][j], aCf[surv & c][j]), c(eC[surv & c][j], eCf[surv & c][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(aBf[surv & c][j], aCf[surv & c][j]), c(eBf[surv & c][j], eCf[surv & c][j]),
	lty=3, col='black', lwd=.5)
	}	# end j in 1:sum(surv[c])
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(surv[c]))	{
	points(aBf[surv & c][j], eBf[surv & c][j], pch=20, col=allcols[1],cex=.6)
	points(aBf[surv & c][j], eBf[surv & c][j], pch=21, col='black',cex=.6)
	}	# end j in 1:sum(surv[c])
# Plot final destinations of C
	for (j in 1:sum(surv[c]))	{
	points(aCf[surv & c], eCf[surv & c], pch=20, col=allcols[4])
	points(aCf[prox & c], eCf[prox & c], pch=21, col=allcols[2])
	}	# end j in 1:sum(surv[c])
#detach(SumAll)
	}	# end k in 1:2

mtext(outer=TRUE,side=1,line= 1.5,cex=1.2,'Semimajor axis (AU)')
mtext(outer=TRUE,side=2,line= 1.3,cex=1.2,'Eccentricity')
mtext(outer=TRUE,side=3,line=0.75,cex=1.3,'Stellar Orbital Parameters')
mtext(outer=TRUE,side=3,line=-.75,cex=1.1,
	'Initial                              Final')
mtext(outer=TRUE,side=4,line=-.20,cex=1.1,
	'Equal masses               Current masses')

dev.off()



