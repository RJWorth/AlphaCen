###############################################################################
#                                 PLOTS
###############################################################################
### Plot distribution of times
print(paste(prefix,'TimeDistribution.pdf',sep=''))
pdf(paste(prefix,'TimeDistribution.pdf',sep=''),width=4,height=8)
par(mfrow=c(2,1))
hB=hist(log10(tB), breaks=tlevs,col=1:ntbins, main='log(tB)')
hC=hist(log10(tC), breaks=tlevs,col=1:ntbins, main='log(tC)')
dev.off()

########################################################################
### Plot initial params, color coded by survival ########################
### a vs a
#pdf(paste(prefix,'aa.pdf',sep=''))
#plot(aB,aC, pch=pindC, col=ind)
#dev.off()

###############################################################################
### Histograms of all parameters
### set up histogram counts:
haBs = hist(aBf[surv], br=br.aB,plot=F)$counts
haBp = hist(aBf[prox], br=br.aB,plot=F)$counts

heBs = hist(eBf[surv],  br=br.e,plot=F)$counts
heBp = hist(eBf[prox],  br=br.e,plot=F)$counts

hiBs = hist(iBf[surv],  br=br.i,plot=F)$counts
hiBp = hist(iBf[prox],  br=br.i,plot=F)$counts

haCs = hist(aCf[surv], br=br.aC,plot=F)$counts
haCp = hist(aCf[prox], br=br.aC,plot=F)$counts

heCs = hist(eCf[surv],  br=br.e,plot=F)$counts
heCp = hist(eCf[prox],  br=br.e,plot=F)$counts

hiCs = hist(iCf[surv],  br=br.i,plot=F)$counts
hiCp = hist(iCf[prox],  br=br.i,plot=F)$counts

pdf(paste(prefix,'HistAllParameters.pdf',sep=''),width=9,height=6)
par(mfrow=c(2,3))
### plot aB
barplot(rbind(haBp,haBs-haBp), space=0, 
	xlab='Semimajor axis of B (AU)',ylab='Counts', main='aBf')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.aB[(0:(n1/2))*2+1],3) )
legend('topleft',legend=c('Surviving simulations','Proxima-like'),
	fill=grey.colors(2)[2:1])

### plot eB
barplot(rbind(heBp,heBs-heBp), space=0, 
	xlab='Eccentricity of B',ylab='Counts', main='eBf')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.e[(0:(n1/2))*2+1],3) )

### plot diB
barplot(rbind(hiBp,hiBs-hiBp), space=0, 
	xlab='Inclination of B (degrees)',ylab='Counts', main='iBf')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.i[(0:(n1/2))*2+1],3) )

### plot daC
barplot(rbind(haCp,haCs-haCp), space=0, 
	xlab='Semimajor axis of C (AU)',ylab='Counts', main='aCf')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.aC[(0:(n1/2))*2+1],3) )

### plot deC
barplot(rbind(heCp,heCs-heCp), space=0, 
	xlab='Eccentricity of C',ylab='Counts', main='eCf')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.e[(0:(n1/2))*2+1],3) )

### plot diC
barplot(rbind(hiCp,hiCs-hiCp), space=0, 
	xlab='Inclination of C',ylab='Counts', main='iCf')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.i[(0:(n1/2))*2+1],3) )

dev.off()

###############################################################################
### Plot changes
### set up histogram counts:
hdaBs = hist(daB[surv], br=br.daB,plot=F)$counts
hdaB1 = hist(daB[prox1], br=br.daB,plot=F)$counts
hdaB2 = hist(daB[prox2], br=br.daB,plot=F)$counts

hdeBs = hist(deB[surv],  br=br.de,plot=F)$counts
hdeB1 = hist(deB[prox1],  br=br.de,plot=F)$counts
hdeB2 = hist(deB[prox2],  br=br.de,plot=F)$counts

hdiBs = hist(diB[surv],  br=br.di,plot=F)$counts
hdiB1 = hist(diB[prox1],  br=br.di,plot=F)$counts
hdiB2 = hist(diB[prox2],  br=br.di,plot=F)$counts

hdaCs = hist(daC[surv], br=br.daC,plot=F)$counts
hdaC1 = hist(daC[prox1], br=br.daC,plot=F)$counts
hdaC2 = hist(daC[prox2], br=br.daC,plot=F)$counts

hdeCs = hist(deC[surv],  br=br.de,plot=F)$counts
hdeC1 = hist(deC[prox1],  br=br.de,plot=F)$counts
hdeC2 = hist(deC[prox2],  br=br.de,plot=F)$counts

hdiCs = hist(diC[surv],  br=br.di,plot=F)$counts
hdiC1 = hist(diC[prox1],  br=br.di,plot=F)$counts
hdiC2 = hist(diC[prox2],  br=br.di,plot=F)$counts

pdf(paste(prefix,'changes.pdf',sep=''),width=9,height=6)
par(mfrow=c(2,3))
### plot daB
barplot(rbind(hdaB2,hdaB1,hdaBs-hdaB1-hdaB2), space=0, 
	xlab='change in aB (AU)',ylab='Counts', main='aBf-aB')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.daB[(0:(n1/2))*2+1],3) )
legend('topleft',legend=c('Surviving simulations','Proxima-like','Big'),
	fill=grey.colors(2)[3:1])

### plot deB
barplot(rbind(hdeB2,hdeB1,hdeBs-hdeB1-hdeB2), space=0, 
	xlab='change in eB',ylab='Counts', main='eBf-eB')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.de[(0:(n1/2))*2+1],3) )

### plot diB
barplot(rbind(hdiB2,hdiB1,hdiBs-hdiB1-hdiB2), space=0, 
	xlab='change in iB',ylab='Counts', main='iBf-iB')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.di[(0:(n1/2))*2+1],3) )

### plot daC
barplot(rbind(hdaC2,hdaC1,hdaCs-hdaC1-hdaC2), space=0, 
	xlab='change in aC (AU)',ylab='Counts', main='aCf-aC')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.daC[(0:(n1/2))*2+1],3) )

### plot deC
barplot(rbind(hdeC2,hdeC1,hdeCs-hdeC1-hdeC2), space=0, 
	xlab='change in eC',ylab='Counts', main='eCf-eC')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.de[(0:(n1/2))*2+1],3) )

### plot diC
barplot(rbind(hdiC2,hdiC1,hdiCs-hdiC1-hdiC2), space=0, 
	xlab='change in iC',ylab='Counts', main='iCf-iC')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.di[(0:(n1/2))*2+1],3) )

dev.off()

###############################################################################
### Plot inclinations
#pdf(paste(prefix,'inc.pdf',sep=''),width=9,height=6)
#par(mfrow=c(2,3))

#hC   = hist( iC,      br=br.i,plot=F)$counts
#hCs  = hist( iC[surv],br=br.i,plot=F)$counts

#hC2  = hist(iCf,      br=br.i,plot=F)$counts
#hC2s = hist(iCf[surv],br=br.i,plot=F)$counts

#hB2  = hist(iBf,      br=br.i,plot=F)$counts
#hB2s = hist(iBf[surv],br=br.i,plot=F)$counts

### iC, with iC of surviving systems colored in
#barplot(rbind(hCs,hC-hCs), space=0, 
#	xlab='i (degrees)',ylab='Counts', main='Initial iC')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iC, with iCf of surviving systems colored in
#barplot(rbind(hC2s,hC2-hC2s), space=0, 
#	xlab='i (degrees)',ylab='Counts', main='Final iC')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iB, with iBf of surviving systems colored in
#barplot(rbind(hB2s,hB2-hB2s), space=0, 
#	xlab='i (degrees)',ylab='Counts', main='Final iB')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iC of surviving systems as % of total iC in that bin
#barplot(100*(hCs/hC), space=0, 
#	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iCf of surviving systems as % of total iCf in that bin
#barplot(100*(hC2s/hC2), space=0, 
#	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iBf of surviving systems as % of total iBf in that bin
#barplot(100*(hB2s/hB2), space=0, 
#	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
#	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

#dev.off()
#######################################################################
### Side-by-side input and output parameters
pdf(paste(prefix,'ae_io.pdf',sep=''),width=8, height=4)
par(mfrow=c(1,2),mar=c(4,4,1,0))

# Initial distribution, a and e for B and C
# make plot borders
plot(.1,.1, pch='.', 
	xlim=c(min(aB,aBf[surv]),max(c(aC,aCf[surv]))), ylim=c(0,1.1),log='x',
	main='',xlab='',ylab='e')
mtext(1,text='a (AU)',line=2.5,cex=1.2)

# plot B where indB=='red'
#points(aB[indB=='red'], eB[indB=='red'], pch='.', col=indB[indB=='red'])
# plot B where indB isn't 'red'
points(aB[indB!='red'], eB[indB!='red'], pch='.', col=indB[indB!='red'])
# plot B from surviving systems (with black borders)
	for (j in 1:sum(surv))	{
	points(aB[surv][j], eB[surv][j], pch=20, col=allcols[1],cex=.6)
	points(aB[surv][j], eB[surv][j], pch=21, col='black',cex=.6)
	}
# plot all initial C parameters
points(aC, eC, pch='.', col=indC)
# plot C from surviving systems bigger
points(aC[surv], eC[surv], pch=20, col=indC[surv])
# add something to proxima-like systems?
# lines connecting B and C from each surviving sim
	for (j in 1:sum(surv))	{
#	points(aB[prox][j], eB[prox][j], pch=21, col=colors()[257],cex=.7)
#	points(aC[prox][j], eC[prox][j], pch=21, col=colors()[257],cex=.9)
	lines(c(aB[surv][j], aC[surv][j]), c(eB[surv][j], eC[surv][j]),
	lty=3, col='black', lwd=.5)
	}

# add legend
	leg4=as.expression(substitute(phantom(1) == 10^t2, 
		list(t2=log10(tmax))))

	leg5=as.expression(substitute(group("[",list(10^t1, 10^t2),")"), 
		list(t1=log10(tcut),t2=log10(tmax))))

	leg6=as.expression(substitute(phantom(1) < 10^t1 , list(t1=log10(tcut))))
legend('topleft', cex=.8, ncol=2, pch=20,
	col=c('white',allcols[1:3],'white',allcols[4:6]),
	legend=c(	expression(t[B] ~ (yrs)),leg4,leg5,leg6,
				expression(t[C] ~ (yrs)),leg4,leg5,leg6))

### Final distribution, a and e for B and C 
par(mar=c(4,3,1,1))
plot(.1,.1, 
	xlim=c(min(aB,aBf[surv]),max(c(aC,aCf[surv]))), ylim=c(0,1.0),log='x',
	main='',xlab='',ylab='')
mtext(1,text='a (AU)',line=2.5,cex=1.2)
# Plot line for orbit with apocenter = 10,000 AU
polygon(c(a1.prx, a2.prx[length(a2.prx):1]), 
		c(e1.prx, e2.prx[length(a2.prx):1]), 
			col='gray',border='gray')
#lines( a1.prx, e1.prx, lty=2)
#lines( a2.prx[length(a2.prx):1], e2.prx[length(a2.prx):1], lty=2)
for (j in 1:sum(surv))	{
# Draw line from initial to final position
lines(c(aB[surv][j], aBf[surv][j]), c(eB[surv][j], eBf[surv][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[surv][j], aCf[surv][j]), c(eC[surv][j], eCf[surv][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(aBf[surv][j], aCf[surv][j]), c(eBf[surv][j], eCf[surv][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(surv))	{
	points(aBf[surv][j], eBf[surv][j], pch=20, col=allcols[1],cex=.6)
	points(aBf[surv][j], eBf[surv][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(aCf[surv], eCf[surv], pch=20, col=indC[surv])

# Outer text
mtext(2,text='e',line=2.5,cex=1.2, outer=TRUE)


dev.off()
########################################################################
# aB in vs. out
#pdf(paste(prefix,'ai_af.pdf',sep=''),width=7.5, height=7.5)
#par(mfrow=c(2,2))

#plot(aB, abs(aBf), pch='.', col=indB, log='xy',#xlim=c(1,max(aBf, na.rm=T)),
#	main='Final vs. initial semimajor axis for B',
#	xlab='Initial a (AU)', ylab='Final a (AU)')
#points(aB[surv], aBf[surv], pch=20, col=indB[surv])
#abline(0,1, lty=2)

#plot(aC, abs(aCf), pch='.', col=indC, log='xy',xlim=c(1,max(aCf, na.rm=T)),
#	main='Final vs. initial semimajor axis for C',
#	xlab='Initial a (AU)', ylab='Final a (AU)')
#points(aC[surv], aCf[surv], pch=20, col=indC[surv])
#abline(0,1, lty=2)

#plot(aB, abs(rBf), pch='.', col=indB, log='xy',#xlim=c(1,max(rBf, na.rm=T)),
#	main='Final position vs. initial semimajor axis for B',
#	xlab='Initial a (AU)', ylab='Final r (AU)')
#points(aB[surv], rBf[surv], pch=20, col=indB[surv])
#abline(0,1, lty=2)

#plot(aC, abs(rCf), pch='.', col=indC, log='xy',xlim=c(1,max(rCf, na.rm=T)),
#	main='Final position vs. initial semimajor axis for C',
#	xlab='Initial a (AU)', ylab='Final r (AU)')
#points(aC[surv], rCf[surv], pch=20, col=indC[surv])
#abline(0,1, lty=2)

#dev.off()
#########################################################################
# Final distribution, a and e for B and C, only where C moved outward
#pdf(paste(prefix,'C_growth.pdf',sep=''),width=4, height=4)
#plot(.1,.1, 
#	xlim=c(1e1,max(c(aC,aCf[grow]))), ylim=c(0,1.0),log='x',
#	main='Output Parameters',xlab='a (AU)',ylab='e')
#for (j in 1:sum(grow))	{
# Draw line from initial to final position
#lines(c(aB[grow][j], aBf[grow][j]), c(eB[grow][j], eBf[grow][j]),
#	col=allcols[1], lwd=.5)
#lines(c(aC[grow][j], aCf[grow][j]), c(eC[grow][j], eCf[grow][j]),
#	col=allcols[5], lwd=.5)
# Connect B and C from each sim
#lines(c(aBf[grow][j], aCf[grow][j]), c(eBf[grow][j], eCf[grow][j]),
#	lty=3, col='black', lwd=.5)
#	}
	# Plot final destinations of B (with fancy borders!)
#	for (j in 1:sum(grow))	{
#	points(aBf[grow][j], eBf[grow][j], pch=20, col=allcols[1],cex=.6)
#	points(aBf[grow][j], eBf[grow][j], pch=21, col='black',cex=.6)
#	}
# Plot final destinations of C
#points(aCf[grow], eCf[grow], pch=20, col=indC[grow])

# make lines showing approximately where C is supposed to be now
#abline(v=)
#abline(h=)

#dev.off()

#########################################################################
# Plot various parameters against time, to see if any tend to live longer
#pdf(paste(prefix,'time.pdf',sep=''), width=8, height=8)
#par(mfrow=c(2,2))

#plot(tB, aB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(aC)), log='x',
#	main='Time vs. initial a', xlab='Time (yrs)',ylab='a (AU)')
#points(tC, aC, pch='.', col=indC)
#points(tB[surv], aB[surv], pch=20, col=indB[surv])
#points(tC[surv], aC[surv], pch=20, col=indC[surv])

#plot(tB, eB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(eC)), log='x',
#	main='Time vs. initial e', xlab='Time (yrs)',ylab='e')
#points(tC, eC, pch='.', col=indC)
#points(tB[surv], eB[surv], pch=20, col=indB[surv])
#points(tC[surv], eC[surv], pch=20, col=indC[surv])

#plot(tB, iB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(iC)), log='x',
#	main='Time vs. initial i', xlab='Time (yrs)',ylab='i (degrees)')
#points(tC, iC, pch='.', col=indC)
#points(tB[surv], iB[surv], pch=20, col=indB[surv])
#points(tC[surv], iC[surv], pch=20, col=indC[surv])

#plot.new()
#legend('bottomright', cex=.7, ncol=2, col=allcols,
#	pch=20,
#	legend=c('B = 1e10','1e5 < B < 1e10','B < 1e5',
#			 'C = 1e10','1e5 < C < 1e10','C < 1e5'))

#dev.off()

###############################################################################

pB=(1-eB)*aB
pC=(1-eC)*aC
pB2=(1-eBf)*aBf
pC2=(1-eCf)*aCf
apB=(1+eB)*aB
apC=(1+eC)*aC
apB2=(1+eBf)*aBf
apC2=(1+eCf)*aCf
peris=cbind(apB,apC,apB2,apC2)
apos =cbind(pB,pC,pB2,pC2)
#inpars=cbind(SumAll[,c(2:3,5:7)],pB,apB,pC,apC)
#outpars=cbind(SumAll[,c(8:13)],pB2,apB2,pC2,apC2)
inpars=cbind(aB,eB,pB,apB,aC,eC,iC,pC,apC)
outpars=cbind(EBf,rBf,aBf,eBf,pB2,apB2,ECf,rCf,aCf,eCf,pC2,apC2)

#panel.time
p.fate=function(x,y)	{points(x,y, pch='.',cex=0.5, col=dcol)}
p.time=function(x,y)	{points(x,y, pch='.',cex=0.5, col=tcol)}

p.fate2=function(x,y)	{points(x,y, pch='.',cex=0.5, col=dcol[surv])}
p.time2=function(x,y)	{points(x,y, pch='.',cex=0.5, col=tcol[surv])}

# pairs
pdf(paste(prefix,'pairs.pdf',sep=''),width=10.5,height=8)
pairs(cbind(aB,eB,aC,eC,iC,EBf,rBf,aBf,eBf,ECf,rCf,aCf,eCf),
	upper.panel=p.fate, lower.panel=p.time)
#dev.off()
pairs(cbind(pB,apB,pC,apC,EBf,aBf,eBf,ECf,aCf,eCf)[surv,],
	upper.panel=p.fate2, lower.panel=p.time2)

#pdf(paste(prefix,'pairs-time.pdf',sep=''),width=10.5,height=8)
#pairs(SumAll[,c(2:3,5:13)],pch='.',cex=0.5,col=tcol)
#dev.off()

#pdf(paste(prefix,'inpairs.pdf',sep=''),width=10.5,height=8)
pairs(inpars,	upper.panel=p.fate, lower.panel=p.time)
#dev.off()

#pdf(paste(prefix,'outpairs.pdf',sep=''),width=10.5,height=8)
pairs(outpars[surv,],	upper.panel=p.fate2, lower.panel=p.time2)
dev.off()
###############################################################################
### Really huge pairs plot!
p.fate2=function(x,y)	{points(x,y, pch=pt$pchs,cex=0.5, col=pt$cols)}

pdf(paste(prefix,'hugepairs.pdf',sep=''),width=24,height=24)
pairs(SumAll[,c(-4,-19,-21:-23)],
	upper.panel=p.fate2, lower.panel=p.time)
dev.off()
###############################################################################

