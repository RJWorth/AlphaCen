# Make all four disk survival plots on one 
###############################################################################

print('Making multiplots')
options(scipen=0)

for (iter in 1:2)	{
if (iter==1)	{
pdf('../Paper/Inserts/DiskSurv.pdf',height=4.5,width=5)
	} else {
png('../Paper/Inserts/DiskSurv.png',height=4.5,width=5,units="in",res=150)
	}
#pdf(paste(simdir,'/DiskSurv.pdf',sep=''), height=7,width=5)
##setEPS()
##postscript(paste('../Paper/Inserts','/DiskSurv.eps',sep=''), height=3.75,width=5)
par(mfrow=c(2,1), oma=c(3,2,0,0) )

marlist=list( c(0.5, 2., 0.5, 1.5), 
	          c(0.5, 2., 0.5, 1.5) )

for (i in 1:length(subdirs))	{

	par(mar=marlist[[i]])

	disk     = l.disk[[i]]  
	diskimg  = l.diskimg[[i]]
	time     = l.time[[i]]
	surv.per = l.surv.per[[i]]
	stab.per = l.stab.per[[i]]
	edge     = l.edge[[i]]

### Plot background image
image(diskimg,col=c(heat,'lightgray','darkgray'), axes=FALSE)
#axis(4, col='red',col.axis="red", lwd=2,las=1,
#	labels = AxLabs,
#	at     = AxLocs)
#mtext(4,text="Disk particle semimajor axis (AU)",line=2.25,col='red')
SimID = paste(substr(    simdir, nchar(    simdir)-1, nchar(    simdir) ),'/',
	      substr(subdirs[i], nchar(subdirs[i])-2, nchar(subdirs[i]) ),sep='' )
#plottitle = paste(SimID,': ',stab.per[length(time)]*100,'% of disk remains stable',sep='')
plottitle = 'Disk Stability: Binary vs. Triple Systems'

### Plot lines over the top
### (xaxs='i' makes axes meet at 0)
par(new=T)
plot(time[-1],r0[edge[-1]]  , lwd=1, lty=1,col='black',
	type='l',log='x',
	xaxs='i',yaxs='i',
	xaxt='n',yaxt='n',
	xlim=c(min(time[-1]),max(time)),ylim=c(0,max(r0)),
	ylab='',xlab='',
	main='')
#lines(time,surv.per, lwd=2)
#lines(time,stab.per, lwd=3, lty=3,col='blue')
#lines(time,edge/n  , lwd=2, lty=2,col='black')

### Add labels to the outer edges of plot
mtext(side=2,'Disk Radius (AU)',       outer=TRUE, line=0.5,cex=1.25,col='black')
#mtext(side=2,'Fraction of disk',       outer=TRUE, line=0.5,cex=1.25,col='black')
#mtext(side=4,expression('Disk particle a/a'[bin]*''),
#	outer=TRUE, line=0.75,cex=1.25,col='red')

### Depending on whether it's the top or bottom plot, add axes/labels
axis(2,lwd=2)
if (i==1)	{
axis(1,labels=FALSE)
#axis(2,at=c(0., .2, .4, .6, .8, 1.0), lwd=2, col='blue',col.axis='blue',las=1)
#mtext(2,text='Surviving fraction', line=2.25)
#mtext(3,text='Disk Stability Without and With Proxima', line=1,cex=1.25)
#mtext(3,text=plottitle, 
#	line=1,cex=1.25)
legend('bottomleft',bty='n',
	legend=c('Unstable','Ejected/accreted'),
	pch=c(15,15),pt.cex=2,
	lty=c(NA,NA),col=c('lightgray','darkgray'),lwd=3)
	} else {
axis(1)
mtext(1,text='Time (yrs)', line=2.25, cex=1.25)
#axis(2,     at=c( 0.0,  0.2,  0.4,  0.6,  0.8, 1.0),
#	labels=c('0.0','0.2','0.4','0.6','0.8', ''), 
#	col='blue',col.axis='blue',lwd=2,las=1)
legend('bottomleft',bty='n',
	legend=c('Truncation Radius'),
	pch=c(NA),pt.cex=2,
	lty=c(1),col=c('black'),lwd=1)
#mtext(2,text='Surviving fraction', line=2.25)

	}
	}	# i in 4

dev.off()
	}	#iter

