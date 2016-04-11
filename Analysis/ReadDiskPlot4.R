# Make all four disk survival plots on one 
###############################################################################

print('Making multiplots')
options(scipen=0)

imgsize = c(4.5,5,150)
#imgname = '../Paper/Inserts/DiskSurv.'
imgname = paste(simdir,'/DiskSurvA.',sep='')
print(imgname)

### Make ps, pdf, and/or png versions of this plot
for (iter in c(3))	{
if (iter==1)	{
setEPS(horizontal=F, onefile=F, paper='special')
postscript(paste(imgname,'eps',sep=''),height=imgsize[1],width=imgsize[2])
	} else if (iter==2) {
pdf(paste(imgname,'pdf',sep=''),height=imgsize[1],width=imgsize[2])
	} else {
png(paste(imgname,'png',sep=''),height=imgsize[1],width=imgsize[2],
									units="in",res=imgsize[3])
	}
### Two plots stacked vertically
par(mfrow=c(2,1), oma=c(3,2,0,0) )

### Top vs bottom plot margins
marlist=list( c(0.5, 2., 0.5, 1.5), 
	          c(0.5, 2., 0.5, 1.5) )

### Loop over subdirectories
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
SimID = paste(substr(    simdir, nchar(    simdir)-1, nchar(    simdir) ),'/',
	      substr(subdirs[i], nchar(subdirs[i])-2, nchar(subdirs[i]) ),sep='' )
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

### Depending on whether it's the top or bottom plot, add axes/labels
axis(2,lwd=2)
if (i==1)	{
axis(1,labels=FALSE)

### Make plot legend
legend('bottomleft',bty='n',
	legend=c('Unstable','Ejected/accreted'),
	pch=c(15,15),pt.cex=2,
	lty=c(NA,NA),col=c('lightgray','darkgray'),lwd=3)
	} else {
axis(1)
mtext(1,text='Time (yrs)', line=2.25, cex=1.25)
legend('bottomleft',bty='n',
	legend=c('Truncation Radius'),
	pch=c(NA),pt.cex=2,
	lty=c(1),col=c('black'),lwd=1)

	}
	}	# i in 4

dev.off()
	}	#iter

