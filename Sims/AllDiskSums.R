
options(width=120)
if (exists('AllSums')) detach(AllSums)

### DiskSummary.txt tags to read in
suf=c('Equal','0717','0723','0814')

### Read in each DiskSummary.txt file and concatenate them
for (i in suf)	{
	if (i=='Equal') m='Equal' else m='True'
	fname=paste('Proxlike/Plots/DiskSummary-',i,'.txt',sep='')
	ThisTab = read.table(fname,header=T)

	Masses=rep(m,dim(ThisTab)[1])
	Dir=rep(i,dim(ThisTab)[1])
	Sim=rownames(ThisTab)
	Reform = cbind( Masses, Dir, Sim, ThisTab )
	print(dim(Reform))

	if (i==suf[1]) AllSums=Reform else AllSums = rbind(AllSums,Reform)
	}
### Format concatenated data array
rownames(AllSums)=1:dim(AllSums)[1]
AllSums$prox=as.factor(AllSums$prox)
AllSums$match=as.factor(AllSums$match)

### Delta(r.a)
AllSums$dR.a = (AllSums$r2.a - AllSums$r3.a)
AllSums$d.frc = AllSums$dR.a/AllSums$r2.a
AllSums$iM   = abs(AllSums$iB - AllSums$iC)
AllSums$pB = AllSums$aB*(1-AllSums$eB)
AllSums$pC = AllSums$aC*(1-AllSums$eC)
#AllSums = cbind(
#	AllSums[,c( -which(colnames(AllSums)=='iB'),
#				-which(colnames(AllSums)=='iC') )],
#	iM, dR.a, d.frc)

rm(Masses)
rm(Dir)
rm(Sim)
#rm(dR.a)
#rm(d.frc)
#rm(iM)

attach(AllSums)


### Index to select for valid simulations only
### better than imported 'prox', which doesn't check time
p = (prox==1) & (t >= 1e7) & (aB > 0)

### Max and min of both deltas
mindR = min( dR.a[p])
maxdR = max( dR.a[p])

minD  = min(delta[p])
maxD  = max(delta[p])

###############################################################################
### Histograms of delta for prox-like True vs. Equal cases
n=15
brD  =  minD + (0:n)*( maxD-minD )/n
brdR = mindR + (0:n)*(maxdR-mindR)/n

pdf('Proxlike/Plots/Hist.pdf')
par(mfrow=c(2,2),oma=c(0,0,1,0))

for (case in c('True','Equal'))	{
hist(delta[prox==1 & Masses==case], 
	br=brD,  col='darkgray', xlab='Added disk truncation (AU)',
	main=paste(case,'masses'))
	}
mtext('Portion of Disk Removed by Proxima',side=3,outer=T,line=-.75,cex=1.5)

for (case in c('True','Equal'))	{
hist( dR.a[prox==1 & Masses==case], 
	br=brdR, col='darkgray', 
	xlab=expression('Added disk truncation (a'[bin]*')'),
	main=paste(case,'masses'))
	}
dev.off()

### Pairs plot
#source('PairsPanels.R')
require(psych)

pairscols = c('match', 'r2.a','r3.a', 'dR.a','d.frc',
	'aB','eB','pB', 'da','de', 'iM', 'aC','eC','pC')
#pairs(AllSums[p,pairscols], pch=20, col=Masses[p])
#pairs.panels(, pchs=20, cols=Masses[p])

pdf('Proxlike/Plots/Pairs.pdf',width=10,height=10)
pairs.panels(AllSums[p,pairscols], scale=T, #jiggle=T,factor=5,
	smooth=F,ellipses=F,lm=F,
	bg=c("black",'red')[AllSums$Mass], hist.col="darkgrey",
	pch=21,
	main="All Proxlike disk sims")
dev.off()

### Get correlations and significances
require(Hmisc)
c=rcorr(as.matrix(AllSums[p,pairscols]))
#print(c$r[c$P>0.05 & c$r>0.5])	# need to do this but keep in matrix form

###
fit1 = lm(r3.a[p] ~ pB[p])
fit2 = lm(r3.a[p]*aB[p] ~ pB[p])

AllSums$fit1.res = rep(NA,dim(AllSums)[1])
AllSums$fit1.res[p] = fit1$residuals
pairscols2 =c('match', 'r2.a','r3.a', 'dR.a','d.frc',
	'aB','eB','pB', 'da','de', 'iM', 'aC','eC','pC','fit1.res')
pairs.panels(AllSums[p,pairscols2], scale=T, #jiggle=T,factor=5,
	smooth=F,ellipses=F,lm=F,
	bg=c("black",'red')[AllSums$Mass], hist.col="darkgrey",
	pch=21,
	main="All Proxlike disk sims")

###############################################################################
### Find medians for each parameter for each mass case
meds = cbind( rep(0,dim(AllSums)[2]),rep(0,dim(AllSums)[2])  )
	colnames(meds) = c('True','Equal')
	rownames(meds)=colnames(AllSums)
for (i in rownames(meds))	{
	for (j in colnames(meds))	{
		meds[i,j] = median(as.numeric(AllSums[p & Masses==j,i]))
	}	}

###############################################################################
### Find 'most median' system(s)
i=4

subset=AllSums[Masses=='Equal' & p==T,]

count = rep(0,dim(subset)[1])
for (i in c(4:8,11:24))	{
mid=which.min( abs(subset[, colnames(subset)[i]] - median( subset[,colnames(subset)[i]]) )) 
print(paste(colnames(subset)[i],mid))
count[mid]=count[mid]+1
	}

###############################################################################
### Compare rtr vs. p for B-2 vs B-3 -- normalized to aBin
pB2=((aB+da)*(1-(eB+de)))[p]
pB3=pB[p]

f2.a=lm(r2.a[p]~pB2)
f3.a=lm(r3.a[p]~pB3)
f2=lm(r2[p]~pB2)
f3=lm(r3[p]~pB3)

for (iter in 1:2)	{
if (iter==1)	{
pdf('../Paper/Inserts/rTra-vs-peri.pdf',height=4.5,width=4.5)
	} else {
png('../Paper/Inserts/rTra-vs-peri.png',height=4.5,width=4.5, units="in",res=150)
	}
plot(pB2,r2.a[p], pch=20,
	xlab='Pericenter (AU)',
	ylab=expression('Truncation radius (r'[tr]*'/a'[bin]*')'),
	xlim=c(min(c( pB2, pB3)),max(c( pB2, pB3))),
	ylim=c(min(c(r2.a[p],r3.a[p])),max(c(r2.a[p],r3.a[p]))) )
abline(f2.a)
points(pB3,r3.a[p], pch=20,col='red')
abline(f3.a,col='red')
legend('topleft',legend=c('Binary','Triple'),
	lty=1,pch=20,col=c('black','red'))
dev.off()
###############################################################################
### Compare rtr vs. p for B-2 vs B-3 -- rTr in AU

if (iter==1)	{
pdf('../Paper/Inserts/rTr-vs-peri.pdf',height=4.5,width=4.5)
	} else {
png('../Paper/Inserts/rTr-vs-peri.png',height=4.5,width=4.5, units="in",res=150)
	}
plot(pB2,r2[p], pch=20,
	xlab='Pericenter (AU)',
	ylab=expression('Truncation radius r'[tr]*' (AU)'),
	xlim=c(min(c( pB2, pB3)),max(c( pB2, pB3))),
	ylim=c(min(c(r2[p],r3[p])),max(c(r2[p],r3[p]))) )
abline(f2)
points(pB3,r3[p], pch=20,col='red')
abline(f3,col='red')
legend('topleft',legend=c('Binary','Triple'),
	lty=1,pch=20,col=c('black','red'))
dev.off()
	} # iter

