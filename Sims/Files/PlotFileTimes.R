
pdf('WallVsSimtime.pdf')
plot(10^logt,dt, pch=20, log='xy', xlab='sim time (yrs)',ylab='walltime (s)',col=machine)
points(10^logsimt, realt[1,], pch=8,col='black',cex=1.75)
points(10^logsimt, realt[2,], pch=8,col='red',cex=1.75)
lines(faket,fakewt.c, col='grey', lty=1, lwd=3)
lines(faket,fakewt.s, col='grey', lty=1, lwd=3)
lines(10^timevalues[t.shrt], realt.dbl[[1]], lty=2, col = "black")
lines(10^timevalues[t.shrt], realt.dbl[[2]], lty=2, col = "red")
lines(10^timevalues[t.long], realt.exp[[1]], lty=1, col = "black")
lines(10^timevalues[t.long], realt.exp[[2]], lty=1, col = "red")
#lines(10^timevalues[t.long], realt.lin[[1]], lty=3, col = "black")
#lines(10^timevalues[t.long], realt.lin[[2]], lty=3, col = "red")

### Show fit models
Cform1=bquote( t[wall] == .(c.f[[1]][[1]][1])^t[sim]^.(c.f[[1]][[1]][2]) )
Sform1=bquote( t[wall] == .(c.f[[2]][[1]][1])^t[sim]^.(c.f[[2]][[1]][2]) )

Cform2=bquote( t[wall] == .(c.f[[1]][[2]][1]) %*%    t[sim]^.(c.f[[1]][[2]][2]) )
Sform2=bquote( t[wall] == .(c.f[[2]][[2]][1]) %*%    t[sim]^.(c.f[[2]][[2]][2]) )

Cform3=bquote( t[wall] == .(c.f[[1]][[3]][1])  +  (.(c.f[[1]][[3]][2])*t[sim]) )
Sform3=bquote( t[wall] == .(c.f[[2]][[3]][1])  +  (.(c.f[[2]][[3]][2])*t[sim]) )
form=expression()	# without this, the other formulas don't display correctly

legend('bottomright', legend=c(Cform1,Sform1,Cform2,Sform2,'prediction',form),
	text.col=c('black','red','black','red','black'),
	col=	 c('black','red','black','red','grey'),
	lty=c(2,2,1,1,1),lwd=c(1,1,1,1,3),pch=c(NA))
### Actual legend
legend('topleft',
	legend=c('run on chloe',
			 'run on shapiro','',
			 'actual times',
			 'average times'),
	col=c('black', 'red', 'white', 'black', 'black'),
	lty=c(      1,     1,      NA,      NA,      NA),
	pch=c(     20,    20,      NA,       20,      8),
	pt.cex=c(   1,     1,       0,       1,    1.75))

dev.off()

