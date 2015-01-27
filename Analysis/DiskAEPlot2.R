i.1=1700
i.2=3300
hcols= c( rep('#FF0000FF',i.1-1), 
		  heat.colors( (i.2-i.1)+1),
		  rep('#FFFFFFFF',length(t)-i.2-1) )
rcols= c( rep('#FF0000FF',i.1-1), 
		  rainbow( (i.2-i.1)+1),
		  rep('#FF0000FF',length(t)-i.2-1) )

#plot(p$PrxCen,pch=20)

pdf(paste(readdir,'/xy.pdf',sep=''))

#plot(x$PrxCen,y$PrxCen,pch=20,cex=.5,col=rcols)
plot(x$AlCenB,y$AlCenB,pch=20,cex=.5,col=rcols)
lines(x$AlCenA,y$AlCenA,col='blue')
#lines(x$AlCenB,y$AlCenB,col='orange')

dev.off()

### Compare to B
pdf(paste(readdir,'/BCvsT.pdf',sep=''),width=6,height=9)
par(mfrow=c(3,1))

plot(t,      Mg,pch=20,cex=.5,col=rcols,log='x')
#plot(t,e$PrxCen,pch=20,cex=.5,col=rcols,log='x')
plot(t,a$AlCenB,pch=20,cex=.5,col=rcols,log='x')
plot(t,e$AlCenB,pch=20,cex=.5,col=rcols,log='x')
dev.off()

### Plot all evolutions
pdf(paste(readdir,'/ae.pdf',sep=''))
par(mfrow=c(1,1))

plot(a$AlCenB, e$AlCenB, pch=20, cex=.5, col='black',log='x',
	xlim=c( 1e-1,max(a[,-3],na.rm=T) ),
	ylim=c( 1e-2, 1 ) )
#points(a[,'PrxCen'], e[,'PrxCen'], pch=20, cex=.5, col='blue')
for (i in 3:length(objs))	{
	points(a[,i], e[,i], pch=20, cex=.5, col=heat.colors(n+25)[i-4])
	}
points(a[nt,], e[nt,], pch=21)
abline(v=aout)
abline(h=etop)

dev.off()

