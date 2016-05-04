### Functions and plotting parameters for ReadElem.R to create PltSystems.pdf

### Solar System planet parameters
SS = array(dim=c(4,6), 
	dimnames=list(c('Mercury','Venus','Earth','Mars'),colnames(elem[[1]])))
SS['Mercury',] = c(  0.387098,   0.205630,   7.005, 0.055, NA, NA)
SS['Venus',]   = c(  0.723332, 0.00677323, 3.39458, 0.815, NA, NA)
SS['Earth',]   = c(1.00000261, 0.01671123, 0.00005,    1., NA, NA)
SS['Mars',]    = c(  1.523679,     0.0935,   1.850, 0.107, NA, NA)

### Make a transparent gray color to use for HZ boxes
diskcol = rgb(190, 190, 190, alpha=100, maxColorValue=255)
HZcol   = rgb(176, 226, 255, alpha=100, maxColorValue=255)

### Text for which rtr and sigma values were used in each sim, to add to plots
SimLabs = vector("expression",length(dirs))
for (i in 1:length(dirs)) { 
#	quote1 = bquote( r[tr] == .(rtr[i])*',' ~ c[sigma] == .(sig[i]) )
	SimLabs[i] = substitute(
        expression( alpha == a. ~ ',' ~ c[sigma] == s.),
        list(a. = alpha[i] , s. = sig[i]))[2]
}

### Index for which objects are stable planets, and planetesimals
plind=list()
pmind=list()
for (i in 1:length(dirs))	{
plind[[i]]=(grepl(paste('P[[:digit:]]{1,}',sep=''), rownames(elem[[i]]) ) & 
			elem[[i]][,'a']>0. & 
            elem[[i]][,'mass'] > PlLim)
pmind[[i]]=(grepl(paste('P[[:digit:]]{1,}',sep=''), rownames(elem[[i]]) ) & 
			elem[[i]][,'a']>0. & 
            elem[[i]][,'mass'] < PlLim)
}

### Which SS plts meet the above mass limit
SSbig = which(SS[,'mass'] >  PlLim)
SSsml = which(SS[,'mass'] <= PlLim)

### max and min semimajor axes of planets
alims = c(1,1)
for (i in 1:length(dirs))	{
	if (min(elem[[i]][plind[[i]],'a']) < alims[1]) alims[1]=min(elem[[i]][plind[[i]],'a'])
	if (max(elem[[i]][plind[[i]],'a']) > alims[2]) alims[2]=max(elem[[i]][plind[[i]],'a'])
	}
alims = c(   0, alims[2]*1.1)
elims = c(-0.1,  ceiling(max(plts$e)*10)/10+0.1)

### text position (x SimLabs, y both, x lost.frac)
textpos = c(0.1,.85*elims[2],alims[2]*.95)


###############################################################################
### Planet systems plot functions

# Create each sim's panel
PlotSim= function(i, side, writeax) {
	# create plot frame
	plot(0.,0., pch=NA, xlab='',ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', 
			xlim=c(alims[1],alims[2]), ylim=elims)
	# HZ and original disk
	rect(xleft=0.693,xright=1.241,  ybottom=-1,    ytop=2,    col=HZcol,   lty=0)
	rect(xleft=0.35, xright=rtr[i], ybottom=-0.01, ytop=0.01, col=diskcol, lty=0)
	### Plot planets (> some mass cutoff) black, smaller objs gray
	points(elem[[i]][plind[[i]],'a'], elem[[i]][plind[[i]],'e'], 
			pch=20, cex=3.*elem[[i]][plind[[i]],'mass']**(1./2.))
	points(elem[[i]][pmind[[i]],'a'], elem[[i]][pmind[[i]],'e'], 
			pch=20, cex=3.*elem[[i]][pmind[[i]],'mass']**(1./2.), col='gray')
	
	# X axis, with labels only in last plot
	Axis(side=1, labels=writeax)
	# Y axis, on left or right depending on which column we're in
	if (side == 'R') {
		axis(side=4, at = c(0, .1, .2, .3, .4, .5), las=2,
                                   labels=c('0','','0.2','','0.4',''),cex.axis=.9)
	} else if (side == 'L') {
		Axis(side=2, at = c(0, .1, .2, .3, .4, .5), las=2,
                                   labels=c('0','','0.2','','0.4',''),cex.axis=.9)
	} else { stopifnot(1==0) }
	text(textpos[1],textpos[2], SimLabs[i], adj = c(0.,1))
	text(textpos[3],textpos[2], adj = c(1,1), 
		labels=paste(sprintf("%.1f", round(lost.frac[i]*100,1)),'% mass loss',sep=''))
}

# Plot SS planets
PlotSS = function() {
	par(mar=c(0,1,0,0))
	plot(0.,0., pch=NA, xlab='',ylab='',, xaxt='n', yaxt='n', xaxs='i', yaxs='i',
	xlim=c(alims[1],alims[2]), ylim=c(-0.05,ceiling(max(plts$e)*10)/10+0.05) )

	rect(xleft=0.951,xright=1.676, ybottom=-1, ytop=2, col=HZcol, lty=0)

	points(SS[SSbig,'a'], SS[SSbig,'e'], pch=20, cex=3.*SS[SSbig,'mass']**(1./2.))
	points(SS[SSsml,'a'], SS[SSsml,'e'], pch=20, cex=3.*SS[SSsml,'mass']**(1./2.),
		col='gray')

	Axis(side=1, labels=T)
	Axis(side=2, at = c(0, .1, .2, .3), labels=c('0','','0.2',''),cex.axis=.9,las=2)
	text(textpos[1], textpos[2],'Solar System', adj = c(0.,1))
}

### Put legend in an empty plot
PlotLegend = function() {
	plot(1, type="n", axes=F, xlab="", ylab="")
	legend('bottom', legend=c('1 earth mass','2 earth masses','3 earth masses'), 
		pch=20, pt.cex=3.*c(1,2,3)**(1./2), bg='white')
}
### Outer margin text
MarginText = function() {
	mtext(side=1,line=3,outer=T,'a (AU)')
	mtext(side=2,line=3,outer=T,'e')
#	maintext = bquote( 'Planetary Systems Formed, ' ~ r[tr] == .(rtr[i]) )
#	mtext(side=3,line=1,outer=T,
#		text=maintext )
}

