
print('hi from R!')

runtype=1	#1 for run normally, 2 or 3 for run on saved files

### Read in Summary.out, normal
if(runtype==1) {prefix='Plots/'
SumAll=read.table('SumAll.out', header=T, na.strings='-')	}
if (runtype==2)	{
prefix='Saved/1Gyr_AllI_NoGas/'
SumAll=read.table(paste(prefix,'SumAll.out',sep=''), header=T,na.strings='-')}
if (runtype==3)	{
prefix='Saved/10Gyr_Tides/'
SumAll=read.table(paste(prefix,'SumAll.out',sep=''), header=T,na.strings='-')}


# To correct a read-in bug that affected 8 sims
print(dim(SumAll))
SumAll=SumAll[SumAll$aB<30 & SumAll$aB>22,]
print(dim(SumAll))

attach(SumAll)

#aB2=abs(aB2)
#aC2=abs(aC2)


#lsummary(SumAll)

#use max time to hopefully get length of sims
#tm=max(c(max(tB),max(tC)))
#ncols=log10(tm)-1
#palette(rainbow(4*ncols))

# make color indices
#indB =   ncols-floor(log10(tB)-2)
#	indB[eB2 >= 1.0] = 2
#indC = 2*ncols+floor(log10(tC)-1)
#	indC[eC2 >= 1.0] = 3*ncols-1
#	indC[indC==-Inf] = 2*ncols+1

# Define cut and max times
tcut=1e5
tmax=max(cbind(tB,tC))
# as strings, with +0 or + removed from sci notation
tcutS=sub('\\+0','',tcut)
tmaxS=sub('\\+0','',tmax)
tcutS=sub('\\+','',tcutS)
tmaxS=sub('\\+','',tmaxS)

### Color indices
allcols=c('red','orange','yellow','darkblue','dodgerblue3','cyan')
indB=replicate(length(aB),allcols[1])
	indB[tB<tmax | eB2>1.0]=allcols[2]
	indB[tB<tcut]=allcols[3]
indC=replicate(length(aC),allcols[4])
	indC[tC<tmax | eC2>1.0]=allcols[5]
	indC[tC<tcut]=allcols[6]


#	indB[tB==tm & eB2 <1.0] = 'red'
#	indB[tB <tm | eB2>=1.0] = 'orange'
#	indB[tB <1000]          = 'yellow'
#	indC[tC==tm & eC2 <1.0] = 'blue'
#	indC[tC <tm | eC2>=1.0] = 'cyan'
#	indC[tC <1000]          = 'green'

ind=tB
	ind[ ((tB==tmax & eB2 <1.0) & (tC==tmax & eC2 <1.0)) ]='green'
	ind[ ((tB <tmax | eB2>=1.0) & (tC==tmax & eC2 <1.0)) ]='orange'
	ind[ ((tB==tmax & eB2 <1.0) & (tC <tmax | eC2>=1.0)) ]='red'
	ind[ ((tB <tmax | eB2>=1.0) & (tC <tmax | eC2>=1.0)) ]='grey'

# symbol indices
pindB = replicate(length(tB),21)
	pindB[(tB == tmax) & (eB2 < 1.0)] = 20
pindC = replicate(length(tC),21)
	pindC[(tC == tmax) & (eC2 < 1.0)] = 20

# survival index
#surv=tB	
	surv = (is.na(aB2)==F & is.na(aC2)==F & eC2<1.0)

# growth index (did C move outward?)
grow = (aC2>aC & eC2<1.0)
	grow[is.na(grow)]=FALSE
prox = (aC2*(1+eC2)>10000 & !is.na(aC2))
print(paste("% of sims with no ejection  =",signif(mean(surv)*100,2)))
print(paste("% of sims where C moves out =",signif(mean(grow)*100,2)))
print(paste("% of Proxima-like sims      =",signif(mean(prox)*100,2)))

# change in B?
bgrowth=abs(aB-aB2)/aB
bgrowth[aB2<0]=NA

# Point size index
PtSz = replicate(length(surv),1)
	PtSz[surv==T] = 10

# Change index -- did parameters change by at least 1%?
aBchange = rbind(aB[surv],aB2[surv],(aB2[surv]-aB[surv])/aB[surv]*100)
aCchange = rbind(aC[surv],aC2[surv],(aC2[surv]-aC[surv])/aC[surv]*100)

eBchange = rbind(eB[surv],eB2[surv],(eB2[surv]-eB[surv])/eB[surv]*100)
eCchange = rbind(eC[surv],eC2[surv],(eC2[surv]-eC[surv])/eC[surv]*100)

#iBchange = rbind(iB[surv],iB2[surv],(iB2[surv]-iB[surv])/iB[surv]*100)
#iCchange = rbind(iC[surv],iC2[surv],(iC2[surv]-iC[surv])/iC[surv]*100)
iBchange = rbind(iB[surv],iB2[surv],(iB2[surv]-iB[surv]))
iCchange = rbind(iC[surv],iC2[surv],(iC2[surv]-iC[surv]))

changes=rbind(round(aBchange[3,]),round(aCchange[3,]), round(eBchange[3,]),
	round(eCchange[3,]), round(iBchange[3,]),round(iCchange[3,]))
row.names(changes)=c('aB','aC','eB','eC','iB*','iC*')
#print(changes)

### destination-based color schemes
dcol=rep('red',length(destC))
	dcol[is.na(destB) & is.na(destC) & eC2>1.0]='magenta'
	dcol[destB=='ejected' & is.na(destC)    ]='green'
	dcol[is.na(destB)     & destC=='ejected']='blue'
	dcol[destB=='ejected' & destC=='ejected']='grey'
	dcol[destC=='Center']='yellow'
	dcol[destC=='AlCenB']='orange'
# Make factor version, with levels sorted from most to least common
dcol2=factor(dcol)
dcol2=factor(dcol2,levels=levels(dcol2)[order(-summary(dcol2))])

### C ejection time color index
mintime=4	#minimum log(t) cutoff
n2=log10(tmax)-mintime+1
br=c(0,10^(mintime:log10(tmax)))
	tcol2=cut(tC,breaks=br,labels=1:n2)
	if(n2==6) {tcol1=cut(tC,breaks=br,labels=c('red','yellow','green','cyan',
						'blue','magenta'))} else
	{tcol1=cut(tC,breaks=br,labels=c('red','yellow','green','cyan',
						'blue','purple','magenta'))}
	tcol2[surv]=='black'

########################################################################
### Plot initial params, color coded by survival ########################
### a vs a
pdf(paste(prefix,'aa.pdf',sep=''))
plot(aB,aC, pch=pindC, col=ind)
dev.off()

########################################################################
# aC in vs. out
pdf(paste(prefix,'aa_C.pdf',sep=''),width=3.75, height=3.75)
plot(aC, abs(aC2), pch='.', col=indC, xlim=c(1,max(aC2, na.rm=T)),log='xy',
	main='Final vs. initial semimajor axis for C',
	xlab='Initial a (AU)', ylab='Final a (AU)')
points(aC[surv], aC2[surv], pch=20, col=indC[surv])
abline(0,1, lty=2)

dev.off()

########################################################################
### Side-by-side input and output parameters
pdf(paste(prefix,'ae_io.pdf',sep=''),width=8, height=4)
par(mfrow=c(1,2),mar=c(4,4,1,0))

# Initial distribution, a and e for B and C
# make plot borders
plot(.1,.1, pch='.', 
	xlim=c(1e1,max(c(aC,aC2[surv]))), ylim=c(0,1.0),log='x',
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
# add something to proxima-like systems
	for (j in 1:sum(prox))	{
#	points(aB[prox][j], eB[prox][j], pch=21, col=colors()[257],cex=.7)
#	points(aC[prox][j], eC[prox][j], pch=21, col=colors()[257],cex=.9)
	lines(c(aB[prox][j], aC[prox][j]), c(eB[prox][j], eC[prox][j]),
	lty=3, col='black', lwd=.5)
	}

# add legend
#	leg1=as.expression(substitute(10^t2 ~ yrs, list(t2=log10(tmax))))
#	leg1=''
	leg4=as.expression(substitute(phantom(1) == 10^t2, 
		list(t2=log10(tmax))))

#	leg2=as.expression(substitute(group("[",list(10^t1, 10^t2),")"), 
#		list(t1=log10(tcut),t2=log10(tmax))))
#	leg2=''
	leg5=as.expression(substitute(group("[",list(10^t1, 10^t2),")"), 
		list(t1=log10(tcut),t2=log10(tmax))))

#	leg3=as.expression(substitute(phantom(0) < 10^t1 , list(t1=log10(tcut))))
#	leg3=''
	leg6=as.expression(substitute(phantom(1) < 10^t1 , list(t1=log10(tcut))))
legend('topleft', cex=.8, ncol=2, pch=20,
	col=c('white',allcols[1:3],'white',allcols[4:6]),
	legend=c(	expression(t[B] ~ (yrs)),leg4,leg5,leg6,
				expression(t[C] ~ (yrs)),leg4,leg5,leg6))

### Final distribution, a and e for B and C 
par(mar=c(4,3,1,1))
plot(.1,.1, 
	xlim=c(1e1,max(c(aC,aC2[surv]))), ylim=c(0,1.0),log='x',
	main='',xlab='',ylab='')
mtext(1,text='a (AU)',line=2.5,cex=1.2)
for (j in 1:sum(surv))	{
# Draw line from initial to final position
lines(c(aB[surv][j], aB2[surv][j]), c(eB[surv][j], eB2[surv][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[surv][j], aC2[surv][j]), c(eC[surv][j], eC2[surv][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(aB2[surv][j], aC2[surv][j]), c(eB2[surv][j], eC2[surv][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(surv))	{
	points(aB2[surv][j], eB2[surv][j], pch=20, col=allcols[1],cex=.6)
	points(aB2[surv][j], eB2[surv][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(aC2[surv], eC2[surv], pch=20, col=indC[surv])

# Outer text
mtext(2,text='e',line=2.5,cex=1.2, outer=TRUE)


dev.off()
#########################################################################
# Final distribution, a and e for B and C, only where C moved outward
pdf(paste(prefix,'C_growth.pdf',sep=''),width=8/2.54, height=8/2.54)
plot(.1,.1, 
	xlim=c(1e1,max(c(aC,aC2[grow]))), ylim=c(0,1.0),log='x',
	main='Output Parameters',xlab='a (AU)',ylab='e')
for (j in 1:sum(grow))	{
# Draw line from initial to final position
lines(c(aB[grow][j], aB2[grow][j]), c(eB[grow][j], eB2[grow][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[grow][j], aC2[grow][j]), c(eC[grow][j], eC2[grow][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(aB2[grow][j], aC2[grow][j]), c(eB2[grow][j], eC2[grow][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(grow))	{
	points(aB2[grow][j], eB2[grow][j], pch=20, col=allcols[1],cex=.6)
	points(aB2[grow][j], eB2[grow][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(aC2[grow], eC2[grow], pch=20, col=indC[grow])

# make lines showing approximately where C is supposed to be now
#abline(v=)
#abline(h=)

dev.off()

#########################################################################
# Plot various parameters against time, to see if any tend to live longer
pdf(paste(prefix,'time.pdf',sep=''), width=8, height=8)
par(mfrow=c(2,2))

plot(tB, aB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(aC)), log='x',
	main='Time vs. initial a', xlab='Time (yrs)',ylab='a (AU)')
points(tC, aC, pch='.', col=indC)
points(tB[surv], aB[surv], pch=20, col=indB[surv])
points(tC[surv], aC[surv], pch=20, col=indC[surv])

plot(tB, eB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(eC)), log='x',
	main='Time vs. initial e', xlab='Time (yrs)',ylab='e')
points(tC, eC, pch='.', col=indC)
points(tB[surv], eB[surv], pch=20, col=indB[surv])
points(tC[surv], eC[surv], pch=20, col=indC[surv])

plot(tB, iB, pch='.', col=indB, xlim=c(1e1,tmax), ylim=c(0, max(iC)), log='x',
	main='Time vs. initial i', xlab='Time (yrs)',ylab='i (degrees)')
points(tC, iC, pch='.', col=indC)
points(tB[surv], iB[surv], pch=20, col=indB[surv])
points(tC[surv], iC[surv], pch=20, col=indC[surv])

plot.new()
legend('bottomright', cex=.7, ncol=2, col=allcols,
#	pch=c(20,replicate(ncols-1,21),20,replicate(ncols-1,21)),
	pch=20,
	legend=c('B = 1e10','1e5 < B < 1e10','B < 1e5',
			 'C = 1e10','1e5 < C < 1e10','C < 1e5'))

dev.off()

###############################################################################
### Plot inclinations
pdf(paste(prefix,'inc.pdf',sep=''),width=9,height=6)
par(mfrow=c(2,3))

n1=8
br1=180*(0:n1)/n1
hC=hist(iC,br=br1,plot=F)$counts
hCs=hist(iC[surv],br=br1,plot=F)$counts

hC2 =hist(iC2,br=br1,plot=F)$counts
hC2s=hist(iC2[surv],br=br1,plot=F)$counts

hB2  =hist(iB2,br=br1,plot=F)$counts
hB2s =hist(iB2[surv],br=br1,plot=F)$counts

# iC, with iC of surviving systems colored in
barplot(rbind(hCs,hC-hCs), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Initial iC')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

# iC, with iC2 of surviving systems colored in
barplot(rbind(hC2s,hC2-hC2s), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Final iC')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

# iB, with iB2 of surviving systems colored in
barplot(rbind(hB2s,hB2-hB2s), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Final iB')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)


# iC of surviving systems as % of total iC in that bin
barplot(100*(hCs/hC), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

# iC2 of surviving systems as % of total iC2 in that bin
barplot(100*(hC2s/hC2), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

# iB2 of surviving systems as % of total iB2 in that bin
barplot(100*(hB2s/hB2), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

dev.off()
###############################################################################

pB=(1-eB)*aB
pC=(1-eC)*aC
pB2=(1-eB2)*aB2
pC2=(1-eC2)*aC2
apB=(1+eB)*aB
apC=(1+eC)*aC
apB2=(1+eB2)*aB2
apC2=(1+eC2)*aC2
peris=cbind(apB,apC,apB2,apC2)
apos =cbind(pB,pC,pB2,pC2)
inpars=cbind(SumAll[,c(2:3,5:7)],pB,apB,pC,apC)
outpars=cbind(SumAll[,c(8:13)],pB2,apB2,pC2,apC2)

#panel.time
p.fate=function(x,y)	{points(x,y, pch='.',cex=0.5, col=dcol)}
p.time=function(x,y)	{points(x,y, pch='.',cex=0.5, col=tcol2)}

palette(rainbow(n2))
# pairs
pdf(paste(prefix,'pairs.pdf',sep=''),width=10.5,height=8)
pairs(SumAll[,c(2:3,5:13)],	upper.panel=p.fate, lower.panel=p.time)
#dev.off()

#pdf(paste(prefix,'pairs-time.pdf',sep=''),width=10.5,height=8)
#pairs(SumAll[,c(2:3,5:13)],pch='.',cex=0.5,col=tcol2)
#dev.off()

#pdf(paste(prefix,'inpairs.pdf',sep=''),width=10.5,height=8)
pairs(inpars,	upper.panel=p.fate, lower.panel=p.time)
#dev.off()

#pdf(paste(prefix,'outpairs.pdf',sep=''),width=10.5,height=8)
pairs(outpars,	upper.panel=p.fate, lower.panel=p.time)
dev.off()

###############################################################################

### Print proxima-like systems to a file (aC2>acutoff)
acutoff=5000	# AU
sink(paste(prefix,'proxlike.txt',sep=''))
options(width=300)
print(SumAll[prox,])
options(width=80)
sink()

###############################################################################
detach(SumAll)

### Calculate average runtime per sim
if (runtype==1) source('CalcRunTime.R')











