#############################################################################
### make set of gif images
for (i in DoAEPlot)	{

	### Make images
	jpeg(paste('img/',dir,'/gifimgs/3xy', sprintf("%04d", i),
		'.jpg',sep=''), 
		height=3, width=3, units='in', res=112)
	par(mfcol=c(1,1), mar=c(3,3,1.5,0.5), oma=c(0,0,0,0))
#	par(bg='black',fg='white',
#		col.axis='white',col.lab='white',col.main='white',col.sub='white')

		pcols=c('yellow','orange','red')
	# Plot box, etc
	plot( 0,0, pch=3, col='black', cex=1,	type='p',
		main=paste(sprintf("%2.1f", (t[i]-t[DoAEPlot[1]])/1e3),' kyr',sep=''),
		xlab='',ylab='', xlim=xlimits, ylim=ylimits)
	mtext('X (AU)', side=1, line=2)
	mtext('Y (AU)', side=2, line=2)
		lines(stars3[[2]]$x[range],stars3[[2]]$y[range],col='red')
#		lines(stars3[[1]]$x[range],stars3[[1]]$y[range],col='orange')
		lines(  (CM2-CM3)[range,1],  (CM2-CM3)[range,2],col='orange')
	# star A
	j=1
	ind=i
	points( matrix(-CM3[i,1]), 				# outer black rings
		    matrix(-CM3[i,2]),
		pch=20, col='black', cex=.5)	
	if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
	points( matrix(-CM3[i,1]),				# inner colored circles
		    matrix(-CM3[i,2]),
		pch=20, col=pcols[j], cex=.3)	
	# B 
	j=2
	points( matrix(stars3[[j-1]]$x[i]), 	# outer black rings
		    matrix(stars3[[j-1]]$y[i]),
		pch=20, col='black', cex=.5)	
	if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
	points( matrix(stars3[[j-1]]$x[ind]),	# inner colored circles
		    matrix(stars3[[j-1]]$y[ind]),
		pch=20, col=pcols[j], cex=.3)
	# C 
	if (mode=='triple')	{
		j=3
		points( matrix(stars3[[j-1]]$x[i]), 	# outer black rings
			    matrix(stars3[[j-1]]$y[i]),
			pch=19, col='black', cex=1.3)	
		if (i<3) ind=i else ind=c(i,i-1,i-2)	# trailing dots and
		points( matrix(stars3[[j-1]]$x[ind]),	# inner colored circles
			    matrix(stars3[[j-1]]$y[ind]),
			pch=19, col=pcols[j], cex=.8)	
		}	
	dev.off()
	}

