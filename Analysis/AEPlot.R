###############################################################################
### Pretty pictures of flight
print('Make a-e plots')

### plot of each coordinate over time
pdf(paste('gif/',dir,'/CoordvsT.pdf',sep=''))
par(mfrow=c(3,1))
plot(     t[DoAEPlot],stars3[[2]]$x[DoAEPlot],type='l',col='red')
	lines(t[DoAEPlot],stars3[[1]]$x[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,1],       type='l',col='blue')
plot(     t[DoAEPlot],stars3[[2]]$y[DoAEPlot],type='l',col='red')
	lines(t[DoAEPlot],stars3[[1]]$y[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,2],       type='l',col='blue')
plot(     t[DoAEPlot],stars3[[2]]$z[DoAEPlot],type='l',col='red')
	lines(t[DoAEPlot],stars3[[1]]$z[DoAEPlot],type='l',col='yellow')
	lines(t[DoAEPlot],-CM3[DoAEPlot,3],       type='l',col='blue')
dev.off()

### Make x-y plot square and centered
xmax = max(-CM3[,1],stars3[[1]]$x,stars3[[2]]$x)
ymax = max(-CM3[,2],stars3[[1]]$y,stars3[[2]]$y)
xmin = min(-CM3[,2],stars3[[1]]$y,stars3[[2]]$y)
ymin = min(-CM3[,2],stars3[[1]]$y,stars3[[2]]$y)
midx = mean(c(xmax,xmin))
midy = mean(c(ymax,ymin))
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

