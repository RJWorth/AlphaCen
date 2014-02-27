
# First time
if(exists('SumAll')==FALSE) {
	library(randomForest)
	source('plot.R')
	attach(SumAll)	}


### a-based color schemes
acol1=log10(aB2)
	acol1[!is.na(acol1)]=cut(acol1[!is.na(acol1)], breaks=n2)
acol2=log10(aC2)
	acol2[!is.na(acol2)]=cut(acol2[!is.na(acol2)], breaks=n2)

### e-based color schemes
ecol1=eB2
	ecol1=cut(ecol1, breaks=n2)
ecol2=eC2
	ecol2=cut(ecol2, 
		breaks=c(min(eC2,na.rm=T)-1,1,max(eC2,na.rm=T)+1))

### index for whether eB==eB2 within X%
x=0.01
eq=rep(FALSE,length(eB))
	eq=eB<(eB2+x) & eB>(eB2-x)

### % of each fate for eq=F,T,NA
eqnums=as.numeric(summary(eq)[2:4])
eqfate=cbind(summary(dcol2[!is.na(eq) & eq==TRUE])/eqnums[2],
		summary(dcol2[!is.na(eq) & eq==FALSE])/eqnums[1],
		summary(dcol2[is.na(eq)])/eqnums[3])
	colnames(eqfate)=c('T','F','NA')
print(eqnums)
print(eqfate)
# make list of colors and meanings in logical order for legend
dcollevs=c('red','magenta','green','blue', 'grey', 'yellow', 'orange')
cat=c('survived','hyperbolic','B ejection','C ejected','both ejected',
	'A-C collision','B-C collision')
###############################################################################
fit=data.frame(tcol1,eB,pC)
form=tcol1~eB*pC
#forest=randomForest(form,data=fit)


######################### Make plots ##########################################
palette(rainbow(n2))
pdf(paste(prefix,'eBf-iBf.pdf',sep=''),width=10.5,height=8)
par(mfrow=c(2,2))

### e vs. i, B final, time colors
plot(pch=20,cex=.5,col=tcol2,eB2,iB2)
legend('topright',pch=20,col=c('white',1:n2),
	legend=c('log(tC)<',mintime-1+(1:n2)))

### e vs. i, B final, fate colors
plot(pch=20,cex=.5,col=dcol,eB2,iB2)
for (i in levels(dcol2)) points(eB2[dcol==i],iB2[dcol==i], pch=20,col=i,cex=.5)
legend('topright',col=dcollevs,pch=20, legend=cat)

#legend('topright',pch=20,col=c('white',1:n), 
#	legend=c('log(tC)<',mintime-1+(1:n)))

### time and fate
plot(tB, tC, log='xy', pch=20,col=dcol,cex=.5)
abline(0,1)
for (i in levels(dcol2)) points(tB[dcol==i],tC[dcol==i], pch=20,col=i,cex=.5)

legend('bottomright',col=dcollevs,pch=20, legend=cat)
dev.off()
###############################################################################

pdf(paste(prefix,'eB_i-f.pdf',sep=''),width=10.5,height=8)
par(mfrow=c(2,2))

### eB in vs. out, time and fate colors
plot(pch=20,cex=.5, col=tcol2[eB2<1], eB[eB2<1], eB2[eB2<1],
	xlab=expression('Initial'~e[B]),ylab=expression('Final'~e[B]))
points(pch=20,eB[prox],eB2[prox])
legend('bottomright',cex=.7, pch=20,col=c('white',1:n2),
	legend=c('log(tC)<',mintime-1+(1:n2)))

plot(pch=20,cex=.5, col=dcol[eB2<1], eB[eB2<1], eB2[eB2<1],
	xlab=expression('Initial'~e[B]),ylab=expression('Final'~e[B]))
	points(pch=20,eB[prox],eB2[prox])
#legend('topright',col=dcollevs,pch=20, legend=cat)

### eB i vs. f, fate colors, off and on line:
plot(pch=20,cex=.5, col=dcol[!eq & eB2<1], eB[!eq & eB2<1], eB2[!eq & eB2<1], 
	main=paste(signif(eqfate[3,2]*100,2),'% survived, ',
				signif(eqfate[1,2]*100,2),'% ejected',sep=''),
	xlab=expression('Initial'~e[B]),ylab=expression('Final'~e[B]))

plot(pch=20,cex=.5, col=dcol[eq], eB[eq], eB2[eq], 
	main=paste(signif(eqfate[3,1]*100,2),'% survived, ',
				signif(eqfate[1,1]*100,2),'% ejected',sep=''),
	xlab=expression('Initial'~e[B]),ylab=expression('Final'~e[B]))
legend('topleft',col=dcollevs,pch=20, legend=cat)

dev.off()
###############################################################################
pdf(paste(prefix,'exploration.pdf',sep=''),width=10.5,height=8)
par(mfrow=c(2,2))


### aB in vs. out
acut=(aB2<50 & aB2>10)
plot(pch=20,cex=.5, col=tcol2[acut], aB[acut], aB2[acut], log='y',
	xlab=expression('Initial'~a[B]),ylab=expression('Final'~a[B]))
lines(c(22,29),c(22,29))
points(pch=20,aB[prox],aB2[prox])
legend('topleft',cex=1, pch=20,col=c('white',1:n2),
	legend=c('log(tC)<',mintime-1+(1:n2)))

#plot(pch=20,cex=.5, col=dcol[acut], aB[acut], aB2[acut], log='y',
#	xlab=expression('Initial'~a[B]),ylab=expression('Final'~a[B]))
#	points(pch=20,aB[prox],aB2[prox])
#legend('topleft',col=dcollevs,pch=20, legend=cat)

plot(pch=20,cex=.5, col=tcol2[eB2<1], eB[eB2<1], eB2[eB2<1],
	xlab=expression('Initial'~e[B]),ylab=expression('Final'~e[B]))
points(pch=20,eB[prox],eB2[prox])


### aC in vs. out
plot(pch=20,cex=.5, col=tcol2[], aC[], aC2[], log='y',
	xlab=expression('Initial'~a[C]),ylab=expression('Final'~a[C]))
points(pch=20,aC[prox],aC2[prox])
legend('topright',cex=1, pch=20,col=c('white',1:n2),
	legend=c('log(tC)<',mintime-1+(1:n2)))

### eC, in vs. out
plot(pch=20,cex=.5, col=tcol2[], eC[], eC2[], log='y',
	xlab=expression('Initial'~e[C]),ylab=expression('Final'~e[C]))
#lines(c(.1,1),c(.1,1))
points(pch=20,eC[prox],eC2[prox])
legend('topright',cex=1, pch=20,col=c('white',1:n2),
	legend=c('log(tC)<',mintime-1+(1:n2)))


### i vs. i, eccentricity colors
#plot(pch=20,cex=.5,col=ecol2,iB2,iC2,main='i from completed sims')
#	points(pch=20,cex=.5,col=ecol2[eC2<1],iB2[eC2<1],iC2[eC2<1])
#legend('topright',pch=20,col=c('white',1:n),legend=c('eC','<1','>1'))

### i vs. i, fate colors (redundant with above?)
#plot(pch=20,cex=.5,col=ind, iB2,iC2)
#legend('topright',pch=20,col=c('grey','red','green'),
#	legend=c('unstable t < 1e8 yrs','unstable in 1e8-1e9','stable until 1e9'))


dev.off()
###############################################################################

par(mfrow=c(1,1))





