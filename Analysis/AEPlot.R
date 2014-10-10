###############################################################################
### Pretty pictures of flight
print('Make a-e plots')

xmax = max(stars[[1]]$x,stars[[2]]$x, na.rm=T)
xmin = min(stars[[1]]$x,stars[[2]]$x, na.rm=T)
ymax = max(stars[[1]]$y,stars[[2]]$y, na.rm=T)
ymin = min(stars[[1]]$y,stars[[2]]$y, na.rm=T)
zmax = max(stars[[1]]$z,stars[[2]]$z, na.rm=T)
zmin = min(stars[[1]]$z,stars[[2]]$z, na.rm=T)
### plot of each coordinate over time
pdf(paste('gif/',dir,'/CoordvsT_Acent.pdf',sep=''))
par(mfrow=c(3,1))
plot(     t[DoAEPlot],stars[[2]]$x[DoAEPlot],type='l',col='red',
															 ylim=c(xmin,xmax))
	lines(t[DoAEPlot],stars[[1]]$x[DoAEPlot],type='l',col='yellow')
#	lines(t[DoAEPlot],-CM3[DoAEPlot,1],       type='l',col='blue')
plot(     t[DoAEPlot],stars[[2]]$y[DoAEPlot],type='l',col='red',
															ylim=c(ymin,ymax))
	lines(t[DoAEPlot],stars[[1]]$y[DoAEPlot],type='l',col='yellow')
#	lines(t[DoAEPlot],-CM3[DoAEPlot,2],       type='l',col='blue')
plot(     t[DoAEPlot],stars[[2]]$z[DoAEPlot],type='l',col='red',
															ylim=c(zmin,zmax))
	lines(t[DoAEPlot],stars[[1]]$z[DoAEPlot],type='l',col='yellow')
#	lines(t[DoAEPlot],-CM3[DoAEPlot,3],       type='l',col='blue')
dev.off()

### Get outermost extent of x and y for all stars
xmax = max(-CM3[,1],stars3[[1]]$x,stars3[[2]]$x, na.rm=T)
xmin = min(-CM3[,1],stars3[[1]]$x,stars3[[2]]$x, na.rm=T)
ymax = max(-CM3[,2],stars3[[1]]$y,stars3[[2]]$y, na.rm=T)
ymin = min(-CM3[,2],stars3[[1]]$y,stars3[[2]]$y, na.rm=T)
zmax = max(-CM3[,3],stars3[[1]]$z,stars3[[2]]$z, na.rm=T)
zmin = min(-CM3[,3],stars3[[1]]$z,stars3[[2]]$z, na.rm=T)

### plot of each coordinate over time
pdf(paste('gif/',dir,'/CoordvsT.pdf',sep=''))
par(mfrow=c(3,1))
plot(     t[DoAEPlot],stars3[[2]]$x[DoAEPlot],type='l',col='red',
															 ylim=c(xmin,xmax))
	lines(t[DoAEPlot],stars3[[1]]$x[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,1],       type='l',col='blue')
plot(     t[DoAEPlot],stars3[[2]]$y[DoAEPlot],type='l',col='red',
															ylim=c(ymin,ymax))
	lines(t[DoAEPlot],stars3[[1]]$y[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,2],       type='l',col='blue')
plot(     t[DoAEPlot],stars3[[2]]$z[DoAEPlot],type='l',col='red',
															ylim=c(zmin,zmax))
	lines(t[DoAEPlot],stars3[[1]]$z[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,3],       type='l',col='blue')
dev.off()

### Make x-y plot square and centered
midx = mean(c(xmax,xmin), na.rm=T)
midy = mean(c(ymax,ymin), na.rm=T)
maxrange=max( xmax-xmin, ymax-ymin )/2
xlimits=c(midx-maxrange,midx+maxrange)
ylimits=c(midy-maxrange,midy+maxrange)
### single line plot of orbit
pdf(paste('gif/',dir,'/XY.pdf',sep=''))
plot(stars3[[2]]$x[DoAEPlot],stars3[[2]]$y[DoAEPlot],type='l',col='red',
	 xlim=xlimits, ylim=ylimits,
	 xlab='X (AU)', ylab='Y (AU)')
	lines(stars3[[1]]$x[DoAEPlot],stars3[[1]]$y[DoAEPlot],col='yellow')
	lines(       -CM3[DoAEPlot,1],       -CM3[DoAEPlot,2],col='blue')
	legend('topleft',c('A','B','C'),col=c('blue','yellow','red'),lty=1)
dev.off()

### Gif images of binary's orbit
source('binary.R')

### Gif images of triple's orbit
if (mode=='triple') source('triple.R')

