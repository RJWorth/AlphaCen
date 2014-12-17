# Make all four disk survival plots on one 
###############################################################################

#pdf(paste(simdir,'/DiskSurv.pdf',sep=''), height=7,width=5)
setEPS()
postscript(paste(simdir,'/DiskSurv.eps',sep=''), height=3.75,width=5)
par(mfrow=c(2,1), oma=c(0,1.5,0,1.5) )

marlist=list( c(0.5, 2., 3.5, 2.), 
	      c(3.5, 2., 0.5, 2.) )

for (i in 1:length(subdirs))	{

	par(mar=marlist[[i]])

	disk     = l.disk[[i]]  
	diskimg  = l.diskimg[[i]]
	time     = l.time[[i]]
	surv.per = l.surv.per[[i]]
	stab.per = l.stab.per[[i]]

### Plot background image
image(diskimg,col=c(heat,'lightgray','darkgray'), axes=FALSE)
axis(4, col='red',col.axis="red", lwd=2,
	labels = AxLabs,
	at     = AxLocs)
#mtext(4,text="Disk particle semimajor axis (AU)",line=2.25,col='red')

plottitle = paste(substr(    simdir, nchar(    simdir)-1, nchar(    simdir) ),'/',
		   substr(subdirs[i], nchar(subdirs[i])-2, nchar(subdirs[i]) ),
		   ': ',stab.per[length(time)]*100,'% of disk remains stable',sep='')

### Plot lines over the top
### (xaxs='i' makes axes meet at 0)
par(new=T)
plot(time[-1],surv.per[-1], type='l',lwd=2,col='blue',log='x',
	xaxs='i',yaxs='i',
	xaxt='n',yaxt='n',
	xlim=c(min(time[-1]),max(time)),ylim=c(0,1),
	ylab='',xlab='',
	main='')
#lines(time,surv.per, lwd=2)
lines(time,stab.per, lwd=3, lty=3,col='blue')

### Add labels to the outer edges of plot
mtext(side=2,'Surviving fraction of disk',       outer=TRUE, line=0.25,cex=1.25,col='blue')
mtext(side=4,'Disk particle semimajor axis (AU)',outer=TRUE, line=0.25,cex=1.25,col='red')

### Depending on whether it's the top or bottom plot, add axes/labels
if (i==1)	{
axis(1,labels=FALSE)
axis(2,at=c(0., .2, .4, .6, .8, 1.0), lwd=2, col='blue',col.axis='blue')
#mtext(2,text='Surviving fraction', line=2.25)
mtext(3,text='Disk Stability Without and With Proxima', line=1,cex=1.25)
legend('bottomleft',bty='n',
	legend=c('Unstable','Ejected/accreted'),
	pch=c(15,15),pt.cex=2,
	lty=c(NA,NA),col=c('lightgray','darkgray'),lwd=3)
	} else {
axis(1)
mtext(1,text='Time (yrs)', line=2.25, cex=1.25)
axis(2,     at=c( 0.0,  0.2,  0.4,  0.6,  0.8, 1.0),
	labels=c('0.0','0.2','0.4','0.6','0.8', ''), 
	col='blue',col.axis='blue',lwd=2)
legend('bottomleft',bty='n',
	legend=c('Fraction surviving','Fraction remaining stable'),
	pch=c(NA,NA),pt.cex=2,
	lty=c(1,3),col=c('blue','blue'),lwd=3)
#mtext(2,text='Surviving fraction', line=2.25)

	}
	}	# i in 4

dev.off()

