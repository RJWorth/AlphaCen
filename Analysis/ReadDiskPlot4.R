# Make all four disk survival plots on one 
###############################################################################

pdf(paste(simdir,'/DiskSurv.pdf',sep=''), height=9,width=6.5)
par(mfrow=c(2,1))

for (i in 1:length(subdirs))	{

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
mtext(4,text="Disk particle semimajor axis (AU)",line=2,col='red')

### Plot lines over the top
par(new=T)
plot(time[-1],surv.per[-1], type='l',lwd=2,col='blue',log='x',
	xaxs='i',yaxs='i',
	xlim=c(min(time[-1]),max(time)),ylim=c(0,1),
	ylab='Surviving percentage',xlab='Time (yrs)',
	main=paste(substr(    simdir, nchar(    simdir)-1, nchar(    simdir) ),'/',
		   substr(subdirs[i], nchar(subdirs[i])-2, nchar(subdirs[i]) ),
		   ': ',stab.per[length(time)]*100,'% of disk remains stable',sep=''))
#lines(time,surv.per, lwd=2)
lines(time,stab.per, lwd=3, lty=3,col='blue')
legend('bottomleft',
	legend=c('Fraction surviving','Fraction remaining stable'),
	lty=c(1,3),col='blue',lwd=3)

	}	# i in 4

dev.off()

