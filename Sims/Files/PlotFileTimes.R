
pdf('WallVsSimtime.pdf')
plot(10^logt,dt, pch=20, log='xy', xlab='sim time (yrs)',ylab='walltime (s)',col=machine)
points(10^logsimt, realt[1,], pch=8,col='black',cex=1.75)
points(10^logsimt, realt[2,], pch=8,col='red',cex=1.75)
lines(10^timevalues, realt.exp1[[1]], lty=2, col = "black")
lines(10^timevalues, realt.exp1[[2]], lty=2, col = "red")
lines(10^timevalues, realt.exp2[[1]], lty=1, col = "black")
lines(10^timevalues, realt.exp2[[2]], lty=1, col = "red")

c=summary(doubleexp.model[[1]])$coefficients[,1]
	c=signif( c(10^(10^(c[1])),c[2]), 2)
s=summary(doubleexp.model[[2]])$coefficients[,1]
	s=signif( c(10^(10^(s[1])),s[2]), 2)

### Show fit models
Cform=bquote( t[wall] == .(c[1]) %*% 10^t[sim]^.(c[2]) )
Sform=bquote( t[wall] == .(s[1]) %*% 10^t[sim]^.(s[2]) )
form=expression()

legend('bottomright', legend=c(Cform,Sform,form),
	text.col=c('black',	'red'),lty=c(0),pch=c(NA))
### Actual legend
legend('topleft',
	legend=c('run on chloe',
			 'run on shapiro','',
			 'actual times',
			 'average times',
			 'walltime ~ (sim t)^n',
			 'walltime ~ m^((sim t)^n)'),
	col=c('black', 'red', 'white', 'black', 'black', 'black', 'black'),
	lty=c(      1,     1,      NA,      NA,      NA,       2,       1),
	pch=c(     20,    20,      NA,       20,      8,      NA,      NA),
	pt.cex=c(   1,     1,       0,       1,    1.75,       1,       1))

dev.off()

