# Make all four disk survival plots on one 
###############################################################################

pdf(paste(simdir,'/DiskSurv.pdf',sep=''), height=12,width=12)
par(mfrow=c(2,2))

for (i in 1:length(subdirs))	{

	disk     = l.disk[[i]]  
	diskimg  = l.diskimg[[i]]
	time     = l.time[[i]]
	surv.per = l.surv.per[[i]]
	stab.per = l.stab.per[[i]]

### Plot background image
image(diskimg,col=c(heat,'lightgray','darkgray'), axes=FALSE)
axis(4, col="red", lwd=2,
	labels=c( min(disk[,1,'r']), 5, 10 ),
	at=c(0,.5,1) )

### Plot lines over the top
par(new=T)
plot(time[-1],surv.per[-1], type='l',lwd=2,col='blue',log='x',
	xaxs='i',yaxs='i',
	xlim=c(min(time[-1]),max(time)),ylim=c(0,1),
	ylab='Surviving percentage',xlab='Time (yrs)',
	main=paste(substr(subdirs[i], nchar(subdirs[i])-2, nchar(subdirs[i])),
		   ': ',stab.per[length(time)]*100,'% of disk remains stable',sep=''))
#lines(time,surv.per, lwd=2)
lines(time,stab.per, lwd=3, lty=3,col='blue')
legend('bottomleft',
	legend=c('Fraction surviving','Fraction remaining stable'),
	lty=c(1,3),col='blue',lwd=3)

	}	# i in 4

dev.off()

