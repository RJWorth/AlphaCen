#############################################################################
### make set of gif images
for (i in DoAEPlot)	{

	### Make images
	jpeg(paste('img/',dir,'/gifimgs/2xy', sprintf("%04d", i),
		'.jpg',sep=''), 
		height=3, width=3, units='in', res=112)
	par(mfcol=c(1,1), mar=c(3,3,1.5,0.5), oma=c(0,0,0,0))
#	par(bg='black',fg='white',
#		col.axis='white',col.lab='white',col.main='white',col.sub='white')

		pcols=c('yellow','orange','red')
	# Plot box, etc
	plot( 0,0, pch=3, col='black', cex=1,	type='p',
		main=paste(sprintf("%2.1f", (t[i]-t[DoAEPlot[1]])/1e3),' kyr',sep=''),
		xlab='',ylab='', 
		xlim=c( min(stars2[[1]]$x),max(stars2[[1]]$x) ),
		ylim=c( min(stars2[[1]]$y),max(stars2[[1]]$y) ))
	mtext('X (AU)', side=1, line=2)
	mtext('Y (AU)', side=2, line=2)
	points(stars2[[2]]$x[range],stars2[[2]]$y[range],pch=20,col='red')
	points(stars2[[1]]$x[range],stars2[[1]]$y[range],pch=20,col='orange')
	points(       -CM2[range,1],       -CM2[range,2],pch=20,col='yellow')
	# star A
	j=1
	ind=i
	points( matrix(-CM2[i,1]), 				# outer black rings
		    matrix(-CM2[i,2]),
		pch=20, col='black', cex=.5)	
#	if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
	points( matrix(-CM2[i,1]),				# inner colored circles
		    matrix(-CM2[i,2]),
		pch=20, col=pcols[j], cex=.3)	
	# B 
	j=2
	points( matrix(stars2[[j-1]]$x[i]), 	# outer black rings
		    matrix(stars2[[j-1]]$y[i]),
		pch=20, col='black', cex=.5)	
#	if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
	points( matrix(stars2[[j-1]]$x[ind]),	# inner colored circles
		    matrix(stars2[[j-1]]$y[ind]),
		pch=20, col=pcols[j], cex=.3)
	# C 
	if (mode=='triple')	{
		j=3
		points( matrix(stars2[[j-1]]$x[i]), 	# outer black rings
			    matrix(stars2[[j-1]]$y[i]),
			pch=19, col='black', cex=1.3)	
		if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
		points( matrix(stars2[[j-1]]$x[ind]),	# inner colored circles
			    matrix(stars2[[j-1]]$y[ind]),
			pch=19, col=pcols[j], cex=.8)	
		}	
	dev.off()
	}

