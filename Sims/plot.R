
cat('hi from R!\n')

###1 for run normally, 2 or 3 for run on saved files
runtype=0	

### Calculate average runtime per sim
if ((runtype==0) | (runtype==1)) source('CalcRunTime.R')

### Read in Summary.out, normal
if(runtype==0) {prefix='Plots/'
	SumAll=read.table('SumAll.out', header=T, na.strings='-')	
	SumAll$tB=10^SumAll$logtB
	SumAll$tC=10^SumAll$logtC
		}
if(runtype==1) {prefix='Plots/'
	SumAll=read.table('SumAll.out', header=T, na.strings='-')	
	SumAll061914=read.table('SumAll061914.out', header=T, na.strings='-')	
	SumAll061714=read.table('SumAll061714.out', header=T, na.strings='-')

	SumAll$tB=10^SumAll$logtB
	SumAll$tC=10^SumAll$logtC
		}
if (runtype==2)	{
	prefix='Saved/1Gyr_AllI_NoGas/'
	SumAll=read.table(paste(prefix,'SumAll.out',sep=''),
		header=T,na.strings='-')}
if (runtype==3)	{
	prefix='Saved/10Gyr_Tides/'
	SumAll=read.table(paste(prefix,'SumAll.out',sep=''),
		header=T,na.strings='-')}


# To correct a read-in bug that affected 8 sims
#print(dim(SumAll))
#SumAll=SumAll[SumAll$aB<30 & SumAll$aB>22,]
#print(dim(SumAll))

attach(SumAll)

#rBf=abs(rBf)
#rCf=abs(rCf)


#lsummary(SumAll)

#use max time to hopefully get length of sims
#tm=max(c(max(tB),max(tC)))
#ncols=log10(tm)-1
#palette(rainbow(4*ncols))

# make color indices
#indB =   ncols-floor(log10(tB)-2)
#	indB[eBf >= 1.0] = 2
#indC = 2*ncols+floor(log10(tC)-1)
#	indC[eCf >= 1.0] = 3*ncols-1
#	indC[indC==-Inf] = 2*ncols+1

# Define cut and max times
tcut=1e5
tmax=max(cbind(tB,tC))
tmin=min(cbind(tB,tC))
#tmin=1e4	#minimum t cutoff
# as strings, with +0 or + removed from sci notation
tcutS=sub('\\+0','',tcut)
tmaxS=sub('\\+0','',tmax)
tcutS=sub('\\+','',tcutS)
tmaxS=sub('\\+','',tmaxS)

### Color indices
allcols=c('red','orange','yellow','darkblue','dodgerblue3','cyan')
indB=replicate(length(aB),allcols[1])
	indB[tB<tmax | EBf>0.]=allcols[2]
	indB[tB<tcut]=allcols[3]
indC=replicate(length(aC),allcols[4])
	indC[tC<tmax | ECf>0.]=allcols[5]
	indC[tC<tcut]=allcols[6]


#	indB[tB==tm & eBf <1.0] = 'red'
#	indB[tB <tm | eBf>=1.0] = 'orange'
#	indB[tB <1000]          = 'yellow'
#	indC[tC==tm & eCf <1.0] = 'blue'
#	indC[tC <tm | eCf>=1.0] = 'cyan'
#	indC[tC <1000]          = 'green'

ind=tB
	ind[ ((tB==tmax & EBf <0.) & (tC==tmax & ECf <0.)) ]='green'
	ind[ ((tB <tmax | EBf>=0.) & (tC==tmax & ECf <0.)) ]='orange'
	ind[ ((tB==tmax & EBf <0.) & (tC <tmax | ECf>=0.)) ]='red'
	ind[ ((tB <tmax | EBf>=0.) & (tC <tmax | ECf>=0.)) ]='grey'

# symbol indices
pindB = replicate(length(tB),21)
	pindB[(tB == tmax) & (EBf <0.)] = 20
pindC = replicate(length(tC),21)
	pindC[(tC == tmax) & (ECf <0.)] = 20

### Indices for the outcomes of simulations, with 2nd set including old sims
# error in reading => na values
brkn = is.na(aBf) | is.na(aCf) | is.na(eBf) | is.na(eCf)
# survival index
surv  = (is.na(destB) & is.na(destC) & eCf <1. & !brkn)
#	 surv[is.na( surv)]=FALSE
# survival index
# growth index (did C move outward?)
grow  = (surv  & rCf>aC)
#	 grow[is.na(grow)]=FALSE
#prox = (rCf*(1+eCf)>10000 & !is.na(rCf))
prox  = (surv  & aCf*(1+eCf)>10000)

# too huge (falsely counted as ejection in Mercury)
huge = ((!is.na(destB) & eBf < 1.) | (!is.na(destC) & eCf < 1.))
	huge[is.na(huge)]=FALSE

### get old version data
	if (runtype==1)	{
	surv2 = (is.na(c(rBf,SumAll061914$rBf,SumAll061714$rBf))==F & 
			 is.na(c(rCf,SumAll061914$rCf,SumAll061714$rCf))==F & 
			 c(ECf,SumAll061914$ECf,SumAll061714$ECf) <0.)
		surv2[is.na(surv2)]=FALSE
	grow2 = (surv2 & 
			c(aCf,SumAll061914$aCf,SumAll061714$rCf)>
										c(aC,SumAll061914$aC,SumAll061714$rC))
		grow2[is.na(grow2)]=FALSE
	prox2 = (surv2 & 
			 c(aCf,SumAll061914$aCf,SumAll061714$rCf)>4000 & 
			 !is.na(c(aCf,SumAll061914$aCf,SumAll061714$rCf)))
		}


### Write totals
i=signif(c(mean(brkn),mean(surv),mean(grow),mean(prox))*100,3)
cat(paste(dim(SumAll)[1],'recent simulations\n'))
cat(paste("broken sims (NAs, etc) = ",i[1],'% (',sum(brkn),')\n',sep=''))
cat(paste("sims with no ejection  = ",i[2],'% (',sum(surv),')\n',sep=''))
cat(paste("sims where C moves out = ",i[3],'% (',sum(grow),')\n',sep=''))
cat(paste("Proxima-like sims      = ",i[4],'% (',sum(prox),')\n',sep=''))

if (runtype==1)	{
	nsims=dim(SumAll)[1]+dim(SumAll061914)[1]+dim(SumAll061714)[1]
	cat(paste(nsims,'completed simulations,',dim(SumAll)[1],'recent\n'))
	cat(paste("% of sims with no ejection  =",signif(mean(surv2)*100,2),'\n'))
	cat(paste("% of sims where C moves out =",signif(mean(grow2)*100,2),'\n'))
	cat(paste("% of Proxima-like sims      =",signif(mean(prox2)*100,2),'\n'))
	}

# change in B?
bgrowth=abs(aB-rBf)/aB
bgrowth[rBf<0]=NA

# Point size index
PtSz = replicate(length(surv),1)
	PtSz[surv==T] = 10

### Amount parameters changed
daB = aBf-aB
deB = eBf-eB
diB = iBf-iB

daC = aCf-aC
deC = eCf-eC
diC = iCf-iC

# Change index -- did parameters change by at least 1%?
aBchange = rbind(aB[surv],rBf[surv],(rBf[surv]-aB[surv])/aB[surv]*100)
aCchange = rbind(aC[surv],rCf[surv],(rCf[surv]-aC[surv])/aC[surv]*100)

eBchange = rbind(eB[surv],eBf[surv],(eBf[surv]-eB[surv])/eB[surv]*100)
eCchange = rbind(eC[surv],eCf[surv],(eCf[surv]-eC[surv])/eC[surv]*100)

##iBchange = rbind(iB[surv],iB2[surv],(iB2[surv]-iB[surv])/iB[surv]*100)
##iCchange = rbind(iC[surv],iC2[surv],(iC2[surv]-iC[surv])/iC[surv]*100)
#iBchange = rbind(iB[surv],iB2[surv],(iB2[surv]-iB[surv]))
#iCchange = rbind(iC[surv],iC2[surv],(iC2[surv]-iC[surv]))

changes=rbind(round(aBchange[3,]),round(aCchange[3,]), round(eBchange[3,]),
	round(eCchange[3,]) ) #, round(iBchange[3,]),round(iCchange[3,]))
row.names(changes)=c('aB','aC','eB','eC') #,'iB*','iC*')
#print(changes)

### destination-based color schemes
dcol=rep('black',length(destC))
	dcol[is.na(destB) & is.na(destC) & ECf>0.]='magenta'	# unstable
	dcol[destB=='ejected' & is.na(destC)    ]='blue'		# B ejected
	dcol[is.na(destB)     & destC=='ejected']='red'			# C ejected
	dcol[destB=='ejected' & destC=='ejected']='grey'		# both ejected
	dcol[destC=='Center']='yellow'							# C hit A
	dcol[destC=='AlCenB']='orange'							# C hit B
# Make factor version, with levels sorted from most to least common
dcol2=factor(dcol)
dcol2=factor(dcol2,levels=levels(dcol2)[order(-summary(dcol2))])

### C ejection time color index
tlevs=floor(log10(tmin)):ceiling(log10(tmax))
ntbins=length(tlevs)		# number of time bins
tcol=factor(floor(log10(tC)), levels=tlevs)
#br=c(0, 10^( log10(tmin):(log10(tmax)+1) ))	# boundaries of time bins
#	tcol2=cut(tC,breaks=br,labels=1:(n2+1),right=FALSE)	# tC binned as numbers
	# tC binned as colors:
#	if(n2==6) {tcol1=factor(tcol2,labels=c('red','yellow','green','cyan',
#						'blue','magenta','black'))} else
#	{tcol1=cut(tC,breaks=br,labels=c('red','yellow','green','cyan',
#						'blue','purple','magenta','black'))}

### Define color palette
#palette(rainbow(ntbins))
#palette(sort(gray.colors(ntbins  ),dec=T))
palette(c(sort(heat.colors(ntbins),dec=T)[2:(ntbins)],'black'))

### Histogram breaks spanning domain of each parameter
# number of sections in histogram
n1=20
# max and min of a for each star, including before and after
aBmin = min(aB[surv],aBf[surv],na.rm=T)
aCmin = min(aC[surv],aCf[surv],na.rm=T)
aBmax = max(aB[surv],aBf[surv],na.rm=T)
aCmax = max(aC[surv],aCf[surv],na.rm=T)
# breaks for a, e, and i, and da, de, and di (the change in each)
br.aB  =          aBmin+(         aBmax-         aBmin)*(0:n1)/n1
br.aC  =          aCmin+(         aCmax-         aCmin)*(0:n1)/n1
br.daB = min(daB[surv])+(max(daB[surv])-min(daB[surv]))*(0:n1)/n1
br.daC = min(daC[surv])+(max(daC[surv])-min(daC[surv]))*(0:n1)/n1
br.e   =  1.*(0:n1)/n1
br.de  =  -1.+2.*(0:n1)/n1
br.i   = 180*(0:n1)/n1
br.di  = -180.+2*180*(0:n1)/n1

###############################################################################
### Plot distribution of times
pdf(paste(prefix,'TimeDistribution.pdf',sep=''),width=4,height=8)
par(mfrow=c(2,1))
hB=hist(log10(tB), breaks=tlevs,col=1:ntbins, main='log(tB)')
hC=hist(log10(tC), breaks=tlevs,col=1:ntbins, main='log(tC)')
dev.off()

########################################################################
### Plot initial params, color coded by survival ########################
### a vs a
pdf(paste(prefix,'aa.pdf',sep=''))
plot(aB,aC, pch=pindC, col=ind)
dev.off()

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
hdaBp = hist(daB[prox], br=br.daB,plot=F)$counts

hdeBs = hist(deB[surv],  br=br.de,plot=F)$counts
hdeBp = hist(deB[prox],  br=br.de,plot=F)$counts

hdiBs = hist(diB[surv],  br=br.di,plot=F)$counts
hdiBp = hist(diB[prox],  br=br.di,plot=F)$counts

hdaCs = hist(daC[surv], br=br.daC,plot=F)$counts
hdaCp = hist(daC[prox], br=br.daC,plot=F)$counts

hdeCs = hist(deC[surv],  br=br.de,plot=F)$counts
hdeCp = hist(deC[prox],  br=br.de,plot=F)$counts

hdiCs = hist(diC[surv],  br=br.di,plot=F)$counts
hdiCp = hist(diC[prox],  br=br.di,plot=F)$counts

pdf(paste(prefix,'changes.pdf',sep=''),width=9,height=6)
par(mfrow=c(2,3))
### plot daB
barplot(rbind(hdaBp,hdaBs-hdaBp), space=0, 
	xlab='change in aB (AU)',ylab='Counts', main='aBf-aB')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.daB[(0:(n1/2))*2+1],3) )
legend('topleft',legend=c('Surviving simulations','Proxima-like'),
	fill=grey.colors(2)[2:1])

### plot deB
barplot(rbind(hdeBp,hdeBs-hdeBp), space=0, 
	xlab='change in eB',ylab='Counts', main='eBf-eB')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.de[(0:(n1/2))*2+1],3) )

### plot diB
barplot(rbind(hdiBp,hdiBs-hdiBp), space=0, 
	xlab='change in iB',ylab='Counts', main='iBf-iB')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.di[(0:(n1/2))*2+1],3) )

### plot daC
barplot(rbind(hdaCp,hdaCs-hdaCp), space=0, 
	xlab='change in aC (AU)',ylab='Counts', main='aCf-aC')
	axis(1,at=(0:(n1/2))*2, lab=signif(br.daC[(0:(n1/2))*2+1],3) )

### plot deC
barplot(rbind(hdeCp,hdeCs-hdeCp), space=0, 
	xlab='change in eC',ylab='Counts', main='eCf-eC')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.de[(0:(n1/2))*2+1],3) )

### plot diC
barplot(rbind(hdiCp,hdiCs-hdiCp), space=0, 
	xlab='change in iC',ylab='Counts', main='iCf-iC')
	axis(1,at=(0:(n1/2))*2, lab=signif(  br.di[(0:(n1/2))*2+1],3) )

dev.off()

###############################################################################
### Plot inclinations
pdf(paste(prefix,'inc.pdf',sep=''),width=9,height=6)
par(mfrow=c(2,3))

hC   = hist( iC,      br=br.i,plot=F)$counts
hCs  = hist( iC[surv],br=br.i,plot=F)$counts

hC2  = hist(iCf,      br=br.i,plot=F)$counts
hC2s = hist(iCf[surv],br=br.i,plot=F)$counts

hB2  = hist(iBf,      br=br.i,plot=F)$counts
hB2s = hist(iBf[surv],br=br.i,plot=F)$counts

### iC, with iC of surviving systems colored in
barplot(rbind(hCs,hC-hCs), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Initial iC')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iC, with iCf of surviving systems colored in
barplot(rbind(hC2s,hC2-hC2s), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Final iC')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iB, with iBf of surviving systems colored in
barplot(rbind(hB2s,hB2-hB2s), space=0, 
	xlab='i (degrees)',ylab='Counts', main='Final iB')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iC of surviving systems as % of total iC in that bin
barplot(100*(hCs/hC), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iCf of surviving systems as % of total iCf in that bin
barplot(100*(hC2s/hC2), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

### iBf of surviving systems as % of total iBf in that bin
barplot(100*(hB2s/hB2), space=0, 
	xlab='i (degrees)',ylab='Percent', main='Surviving Systems')
	axis(1,at=n1*(0:4)/4, lab=180*(0:4)/4)

dev.off()
#######################################################################
### Side-by-side input and output parameters
pdf(paste(prefix,'ae_io.pdf',sep=''),width=8, height=4)
par(mfrow=c(1,2),mar=c(4,4,1,0))

# Initial distribution, a and e for B and C
# make plot borders
plot(.1,.1, pch='.', 
	xlim=c(1e1,max(c(aC,rCf[surv]))), ylim=c(0,1.0),log='x',
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
legend('bottomright', cex=.8, ncol=2, pch=20,
	col=c('white',allcols[1:3],'white',allcols[4:6]),
	legend=c(	expression(t[B] ~ (yrs)),leg4,leg5,leg6,
				expression(t[C] ~ (yrs)),leg4,leg5,leg6))

### Final distribution, a and e for B and C 
par(mar=c(4,3,1,1))
plot(.1,.1, 
	xlim=c(1e1,max(c(aC,rCf[surv]))), ylim=c(0,1.0),log='x',
	main='',xlab='',ylab='')
mtext(1,text='a (AU)',line=2.5,cex=1.2)
for (j in 1:sum(surv))	{
# Draw line from initial to final position
lines(c(aB[surv][j], rBf[surv][j]), c(eB[surv][j], eBf[surv][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[surv][j], rCf[surv][j]), c(eC[surv][j], eCf[surv][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(rBf[surv][j], rCf[surv][j]), c(eBf[surv][j], eCf[surv][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(surv))	{
	points(rBf[surv][j], eBf[surv][j], pch=20, col=allcols[1],cex=.6)
	points(rBf[surv][j], eBf[surv][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(rCf[surv], eCf[surv], pch=20, col=indC[surv])

# Outer text
mtext(2,text='e',line=2.5,cex=1.2, outer=TRUE)


dev.off()
########################################################################
# aB in vs. out
pdf(paste(prefix,'ai_af.pdf',sep=''),width=7.5, height=7.5)
par(mfrow=c(2,2))

plot(aB, abs(aBf), pch='.', col=indB, log='xy',#xlim=c(1,max(aBf, na.rm=T)),
	main='Final vs. initial semimajor axis for B',
	xlab='Initial a (AU)', ylab='Final a (AU)')
points(aB[surv], aBf[surv], pch=20, col=indB[surv])
abline(0,1, lty=2)

plot(aC, abs(aCf), pch='.', col=indC, log='xy',xlim=c(1,max(aCf, na.rm=T)),
	main='Final vs. initial semimajor axis for C',
	xlab='Initial a (AU)', ylab='Final a (AU)')
points(aC[surv], aCf[surv], pch=20, col=indC[surv])
abline(0,1, lty=2)

plot(aB, abs(rBf), pch='.', col=indB, log='xy',#xlim=c(1,max(rBf, na.rm=T)),
	main='Final position vs. initial semimajor axis for B',
	xlab='Initial a (AU)', ylab='Final r (AU)')
points(aB[surv], rBf[surv], pch=20, col=indB[surv])
abline(0,1, lty=2)

plot(aC, abs(rCf), pch='.', col=indC, log='xy',xlim=c(1,max(rCf, na.rm=T)),
	main='Final position vs. initial semimajor axis for C',
	xlab='Initial a (AU)', ylab='Final r (AU)')
points(aC[surv], rCf[surv], pch=20, col=indC[surv])
abline(0,1, lty=2)

dev.off()
#########################################################################
# Final distribution, a and e for B and C, only where C moved outward
pdf(paste(prefix,'C_growth.pdf',sep=''),width=8/2.54, height=8/2.54)
plot(.1,.1, 
	xlim=c(1e1,max(c(aC,rCf[grow]))), ylim=c(0,1.0),log='x',
	main='Output Parameters',xlab='a (AU)',ylab='e')
for (j in 1:sum(grow))	{
# Draw line from initial to final position
lines(c(aB[grow][j], rBf[grow][j]), c(eB[grow][j], eBf[grow][j]),
	col=allcols[1], lwd=.5)
lines(c(aC[grow][j], rCf[grow][j]), c(eC[grow][j], eCf[grow][j]),
	col=allcols[5], lwd=.5)
# Connect B and C from each sim
lines(c(rBf[grow][j], rCf[grow][j]), c(eBf[grow][j], eCf[grow][j]),
	lty=3, col='black', lwd=.5)
	}
	# Plot final destinations of B (with fancy borders!)
	for (j in 1:sum(grow))	{
	points(rBf[grow][j], eBf[grow][j], pch=20, col=allcols[1],cex=.6)
	points(rBf[grow][j], eBf[grow][j], pch=21, col='black',cex=.6)
	}
# Plot final destinations of C
points(rCf[grow], eCf[grow], pch=20, col=indC[grow])

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

### Print proxima-like systems to a file (rCf>acutoff)
acutoff=5000	# AU
sink(paste(prefix,'proxlike.txt',sep=''))
options(width=300)
print(SumAll[prox,],row.names=FALSE)
options(width=80)
sink()

###############################################################################
detach(SumAll)











