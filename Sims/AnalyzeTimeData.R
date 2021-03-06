
### Read version read in for first submitted paper draft
fname = 'TimeData/DiskSummary-v1.txt'

hdr = read.table(fname,sep='|',strip.white=T,skip=1,nrows=1, stringsAsFactors=F)
nlines = length(readLines(fname))

disksum = read.table(fname, sep='|', strip.white=T, skip=3, nrows=nlines-4,
	col.names=hdr)

disksum = disksum[,c(-1,-dim(disksum)[2])]
attach(disksum)
# dr vrs tf aBf eBf pBf apBf aCf eCf pCf apCf minpB mintB minpC mintC iMf rtrf rpf rpm

nrows = dim(disksum)[1]
ncols = dim(disksum)[2]
dirnames = dr[1:nrows %% 3 == 0]
stable = (tf>=1e6) & !( (aBf<0 | is.na(aBf)) | (aCf<0 & !is.na(aCf)) )
prx = stable & disksum$apCf>10000 & !is.na(disksum$apCf)
v2 = (vrs=='2')
v3 = (vrs=='3')
Dcase = grepl('d',dr) 

# save list of proxlike sims
listprox = data.frame( prx[vrs=='O'], prx[vrs=='3'],row.names=dirnames)
colnames(listprox) = c('O','B3')
sink('ProxlikeDirs.txt')
cat('# Disk sim directories which satisfy "Proxima-like" requirements\n')
cat('# Made by AnalyzeTimeData.R\n')
print(listprox)
sink()

# Indices for good B-3 dirs:
i3 = v3 & prx & stable
dr3 = dr[i3]
# Corresponding B-2 dirs:
i2 = dr %in% dr3 & v2 & stable
dr2 = dr[i2 & stable]
# circle back to remove B-3s whose B-2s weren't stable
i3 = dr %in% dr2 & v3
dr3 = dr[i3]

### 5 points with r/p ratio outside 0.1-0.4 
outliers = (rtrf/minpB <0.1 | rtrf/minpB>0.4) & !is.na(rtrf) & stable

fitF   = lm(rtrf[stable] ~ pBf[stable])
fitMin = lm(rtrf ~ minpB)
fitcalm= lm(rtrf[v2] ~ minpB[v2])

### PLOTS
#pairs(disksum[,3:ncols],pch=20,col=prx+1)

# good fit
#plot(pBf[stable],rtrf[stable], pch=20, col=prx[stable]+1)
#abline(fitF)
# better
#plot(minpB,rtrf, pch=20, col=prx+1)
#abline(fitMin)

# How much does proxima's interaction change the system over time?
#plot(pBf[stable], minpB[stable], pch=20, col=prx[stable]+1 )
#plot(pBf[prx], minpB[prx], pch=20, col=as.numeric(vrs[prx])-1 )
#plot(pBf[stable], minpB[stable], pch=20-prx[stable], col=as.numeric(vrs[stable])-1 )

#plot(minpB,rtrf,pch=20,col='blue')
#abline(fitMin,col='blue')
#### rtr vs minpB in B-2 systems which were not perturbed:
#points(minpB[v2],rtrf[v2], pch=20)
#abline(fitcalm)
##abline(fitF,col='red')
## the two hig-ratio outlier points pull the 
#fitPrx = lm(rtrf[prx & !outlie] ~ minpB[prx & !outlie])
#points(minpB[prx & !outlie],rtrf[prx & !outlie], pch=20,col='red')
#abline(fitPrx,col='red')

### How much lower is minpB than pBf?
dp = (pBf - minpB)
fp = dp/pBf
PmPf = minpB/pBf

#pericenter reduction from B-2 to B-3
delp  =   pBf[i2 & Dcase] -   pBf[i3 & Dcase]
delpm = minpB[i2 & Dcase] - minpB[i3 & Dcase]
# fractional decrease
fpm   = delpm/minpB[i2 & Dcase]

delpA  =   pBf[i2 & !Dcase] -   pBf[i3 & !Dcase]
delpmA = minpB[i2 & !Dcase] - minpB[i3 & !Dcase]
# fractional decrease
fpmA   = delpmA/minpB[i2 & !Dcase]
p3p2 = minpB[i3]/minpB[i2]
p3p2A = minpB[i3 & !Dcase]/minpB[i2 & !Dcase]
p3p2D = minpB[i3 & Dcase]/minpB[i2 & Dcase]

PfPi = pBf[i3 & Dcase]/pBf[i2 & Dcase]


###############################################################################

### Stats array
vecs = list( rpf2 = rpf[i2], rpf3 = rpf[i3], 
			 rpm2 = rpm[i2], rpm3 = rpm[i3 & is.finite(rpm)], 
			 fp2  =  fp[i2], fp3  =  fp[i3] )
Dvec = list( rpf2 = rpf[i2 & Dcase], rpf3 = rpf[i3 & Dcase], 
			 rpm2 = rpm[i2 & Dcase], rpm3 = rpm[i3 & Dcase & is.finite(rpm)], 
			 fp2  =  fp[i2 & Dcase], fp3  =  fp[i3 & Dcase] )
Avec = list( rpf2 = rpf[i2 & !Dcase], rpf3 = rpf[i3 & !Dcase], 
			 rpm2 = rpm[i2 & !Dcase], rpm3 = rpm[i3 & !Dcase & is.finite(rpm)], 
			 fp2  =  fp[i2 & !Dcase], fp3  =  fp[i3 & !Dcase] )
GetStats = function( vecs) {
	mn  = rep(NA, length(vecs))
	std = rep(NA, length(vecs))
	med = rep(NA, length(vecs))
	mad = rep(NA, length(vecs))

	for (j in 1:length(vecs)) {
		mn[j]  =   mean(vecs[[j]])
		std[j] =     sd(vecs[[j]])
		med[j] = median(vecs[[j]])
		mad[j] =    mad(vecs[[j]])
	}
	stats = t(data.frame(mn,std,med,mad, row.names=names(vecs)))
	return(stats)
}
stats  = GetStats(vecs)
Dstats = GetStats(Dvec)
Astats = GetStats(Avec)

###############################################################################

### 
mSun        = 1.9891e30		# kg
mEarth      = 5.972e24 		# kg
BinPlts = read.table('../PlanetSims/BinPlts.out',header=T)
BinPlts$mass = BinPlts$mass*mSun/mEarth
BinPlts$Sim = as.factor(BinPlts$Sim)
BinPlts = BinPlts[BinPlts$Obj != 'AlCenA',]
attach(BinPlts)
lvls = levels(Sim)
# Add columns for sigma, rtr?
sigmalvls = rep( c(1, 1, 3, 3, .3, .3), 2)
rtrlvls   = rep( c(2.77, 3.08), 6)
sigma = rep(NA, dim(BinPlts)[1])
rtr   = rep(NA, dim(BinPlts)[1])
for (j in 1:dim(BinPlts)[1]) {
	sigma[j] = sigmalvls[ lvls == Sim[j] ]
	rtr[j]   =   rtrlvls[ lvls == Sim[j] ]
}
### Calculate stats about planet sims
nPlts = summary(Sim)-1
HZbounds = c(0.693, 1.241) # AU
HZ = a>=HZbounds[1] & a<=HZbounds[2]
nHZ      = rep(NA, length(nPlts))
MeanMass = rep(NA, length(nPlts))
StDMass  = rep(NA, length(nPlts))
MedMass  = rep(NA, length(nPlts))
MadMass  = rep(NA, length(nPlts))
for (j in 1:length(lvls)) {
	MeanMass[j] =   mean(mass[Sim == lvls[j]])
	StDMass[j]  =     sd(mass[Sim == lvls[j]])
	MedMass[j]  = median(mass[Sim == lvls[j]])
	MadMass[j]  =    mad(mass[Sim == lvls[j]])
	nHZ[j]      =    sum(  HZ[Sim == lvls[j]])
}

ACPStats = data.frame( sigmalvls, rtrlvls,nPlts, nHZ, MedMass,MadMass, MeanMass,StDMass )

#SigStats = data.frame( nPl=c(  ), 
#	nHZ=, 
#	mass=)

###############################################################################

#png('../Paper/Inserts/Rtr_vs_minP.png',height=6,width=6,units="in",res=150)
pdf('../Paper/Inserts/Rtr_vs_minP.pdf',width=5, height=10)
#setEPS(horizontal=F, onefile=F, paper='special')
#postscript('../Paper/Inserts/Rtr_vs_minP.eps',width=5, height=9)

par(mfrow=c(2,1))
par(oma=c( 4.75,    0, 2.5, 1.5))
par(mar=c(    0,  4.5,    0,  0))
# Plot of rtr vs minpB, col = B-2/3, size = Prx or not
plot(minpB[stable], rtrf[stable], 
	pch=20-prx[stable], col=as.numeric(vrs[stable])-1, 
	ylab='',xlab='',main = '')
mtext(side=2,line=2.5,cex=1.5,expression(r[tr]~'(AU)')) 
#fit1  = lm( rtrf[stable] ~ minpB[stable])
#fit2  = lm( rtrf[stable & !outliers] ~ minpB[stable & !outliers])
#fitr1 = lm( rtrf[stable & v3] ~ minpB[stable & v3])
x1 = minpB[stable & v2]
y1 = rtrf[stable & v2]
fitb  = lm( y1 ~ x1)
x2 = minpB[stable & v3 & !outliers]
y2 = rtrf[stable & v3 & !outliers] 
fitr2 = lm(y2 ~ x2)
EqnLabs = vector("expression",2)
	EqnLabs[1] = substitute(expression( r[tr] == a. + b. * p[min]),
		list(a. = round(fitb$coefficients[1],2),
			 b. = round(fitb$coefficients[2],2)) )[2]
	EqnLabs[2] = substitute(expression(r[tr] == a. + b. * p[min]),
		list(a. = round(fitr2$coefficients[1],2),
			 b. = round(fitr2$coefficients[2],2)) )[2]
#abline(fit1, lty=3)
#abline(fit2, lty=2)
abline(fitb, col='black')
text(13,2,adj=c(0,.5),col='black', labels=EqnLabs[1],cex=1.25)
#abline(fitr1, col='red',lty=3)
abline(fitr2, col='red',lty=1)
text(2,6,adj=c(0,.5),col='red', labels=EqnLabs[2],cex=1.25)

plot(minpB[stable], rtrf[stable]/minpB[stable], 
	pch=20-prx[stable], col=as.numeric(vrs[stable])-1, log='y',
	ylab='',xlab='',main = '')
mtext(side=2,line=2.5,cex=1.5,expression(r[tr]/p[bin])) 
mtext(side=1,outer=T,line=3, expression(p[bin]),cex=1.5,adj=c(.6))


dev.off()
###############################################################################

#png('../Paper/Inserts/Rtr_vs_minP.png',height=6,width=6,units="in",res=150)
pdf('../Paper/Inserts/PeriHist.pdf',width=5, height=9)
#setEPS(horizontal=F, onefile=F, paper='special')
#postscript('../Paper/Inserts/PeriHist.eps',width=5, height=9)
par(mfrow=c(2,1), oma=c(0,0,1,0),mar=c(4.5,4.5,0,1))
### Histogram of minimum pericenter as fraction of final pericenter, 
### for 3-star simulations which result in proxima-like systems
hist(PmPf[i3 & Dcase], probability=F, col='gray30', br=10, 
	xlab = expression(p[min]/p[f]), main='',ylab='# of Systems',cex.lab=1.5)
### Reduction in pericenter of 3-star systems which result in Proxiima-like 
### systems, relative to 2-star equivalent, i.e. relative to initial params
hist(p3p2D, probability=F, col='gray30', br=10, 
	xlab = expression(p[min]/p[i]), main='',ylab='# of Systems',cex.lab=1.5)
### Final pericenter as fraction of initial pericenter
#hist(PfPi, probability=T, col='grey', br=10, 
#	xlab = expression(p[f]/p[i]), main='')

dev.off()

###############################################################################
print(paste('Disk Sims: w/equal mass w/o 3rd star, median rtr:',  Dstats['med','rpm2'],'+-',Dstats['mad','rpm2']))
print(paste('Disk Sims: w/equal mass w/  3rd star, median rtr:',  Dstats['med','rpm3'],'+-',Dstats['mad','rpm3']))

print(paste('Min(pB) B-3 smaller than in B-2 by (in %):',median(fpm),'+-',mad(fpm) ))
print(paste('Min(pB) smaller than final pB by (in %):',median(fp[i3 & Dcase]),'+-',mad(fp[i3 & Dcase])))

print(paste('Planets per sim (mn/med):',mean(nPlts),median(nPlts)))
print(paste('HZ Planets per sim (mn/med):',mean(nHZ),median(nHZ)))
print(paste('Avg planet mass, mEarth (mn/med):',mean(MeanMass),mean(MedMass)))

###############################################################################
detach(disksum)
detach(BinPlts)


