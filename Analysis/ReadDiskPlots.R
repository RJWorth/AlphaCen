pdf(paste(dir,'/DiskPairs.pdf',sep=''), height=10,width=10)
pairs(pairframe,pch=20,col=grays)
dev.off()
###
pdf(paste(dir,'/DiskSurv.pdf',sep=''), height=4,width=8)
par(mar=c(4, 4, 4, 4) + 0.1)

image(diskimg,col=c(heat,'lightgray','darkgray'), axes=FALSE)
axis(4, col='red',col.axis="red", lwd=2,
	labels = AxLabs,
	at     = AxLocs)
mtext(4,text="Disk particle semimajor axis (AU)",line=2,col='red')


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
	legend=c('Unstable','Ejected/accreted',
		 'Fraction surviving','Fraction remaining stable'),
	fill=c('lightgray','darkgray','white','white'), border='white',
	lty=c(NA,NA,1,3),col='blue',lwd=3)

dev.off()
###
pdf(paste(dir,'/DiskEvol.pdf',sep=''), height=15,width=15)
par(mfrow=c(3,1))

for (column in c('x','y','z'))	{
plot(time[-1],star[2,time.all[-1],column], type='l', col='orange',log='x',
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(time[-1],disk[i,time.all[-1],column], col=grays[i])
	}

dev.off()
###
pdf(paste(dir,'/DiskOrbits.pdf',sep=''), height=10,width=10)

plot(0,0, pch=19, col='orange', main='End of sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      points(disk[i,timeslice,4:5], col=grays[i],pch=20)
for (i in 2:nstars) points(star[i,timeslice,4:5], col='orange',pch=20)

plot(0,0, pch=19, col='orange', main='Whole sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(disk[i, time.all,4:5], col=grays[i])
for (i in 2:nstars) lines(star[i, time.all,4:5], col='orange')

plot(0,0, pch=19, col='orange', main='Start of sim',
	xlim=c(-extent,extent),
	ylim=c(-extent,extent)	)
for (i in n:1)      lines(disk[i,     1:50,4:5], col=grays[i])
for (i in 2:nstars) lines(star[i,     1:50,4:5], col='orange')
dev.off()

