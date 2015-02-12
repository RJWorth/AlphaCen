
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
### Find 'most median' system(s)
i=4

subset=AllSums[Masses=='Equal' & p==T,]


count = rep(0,dim(subset)[1])
for (i in c(4:8,11:24))	{
mid=which.min( abs(subset[, colnames(subset)[i]] - median( subset[,colnames(subset)[i]]) )) 
print(paste(colnames(subset)[i],mid))
count[mid]=count[mid]+1
	}




