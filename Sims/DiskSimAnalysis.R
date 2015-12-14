suf=c('Equal','0717','0723','0814')

### Read in each DiskSummary.txt file and concatenate them
for (i in suf)	{
	if (i=='Equal') m='Equal' else m='True'
	fname=paste('Proxlike/Plots/DiskSummary-',i,'.txt',sep='')
	ThisTab = read.table(fname,header=T)

	# Make columns with set-identifying information, then combine columns
	Masses=rep(m,dim(ThisTab)[1])
	Dir=rep(i,dim(ThisTab)[1])
	Sim=rownames(ThisTab)

	Reform = cbind( Masses, Dir, Sim, ThisTab )
	print(dim(Reform))

# alternate restructuring attempt that doesn't work yet:
	attach(ThisTab)
	aBi = aB+da
	aBf = aB
	eBi = eB+de
	eBf = eB
	iM  = abs(iB - iC)
	pBi = aBi*(1-eBi)
	pBf = aB*(1-eB)
	pCf = aC*(1-eC)
	detach(ThisTab)
	Reform1 = data.frame(cbind( Masses, Dir, Sim, 
		aBi, eBi, aBf, eBf, pBi, pBf, pCf, iM, ThisTab ))
	remove(Dir)
	remove(Masses)
	remove(Sim)
	remove(aBi)
	remove(aBf)
	remove(eBi)
	remove(eBf)
	remove(iM)
	remove(pBi)
	remove(pBf)
	remove(pCf)
#	# r2, r3, prox, t, aC, eC needed from ThisTab, but not using the whole
#	# thing makes the resulting array all factors, not numerical
	Reform2 = Reform1[,c('Masses','Dir','Sim','r2','r3','prox','t',
	'aBi','eBi','aBf','eBf','aC','eC','iM','pBi','pBf','pCf')]

	# if this isn't the first rows of the table, add it to the previous rows
	if (i==suf[1]) AllSums=Reform2 else AllSums = rbind(AllSums,Reform2)
	}
### Format concatenated data array
rownames(AllSums)=1:dim(AllSums)[1]
AllSums$prox=as.factor(AllSums$prox)
#AllSums$match=as.factor(AllSums$match)

### Delta(r.a)
#AllSums$dR.a = (AllSums$r2.a - AllSums$r3.a)
#AllSums$d.frc = AllSums$dR.a/AllSums$r2.a
#AllSums$iM   = abs(AllSums$iB - AllSums$iC)
#AllSums$pB = AllSums$aB*(1-AllSums$eB)
#AllSums$pC = AllSums$aC*(1-AllSums$eC)

### Index to select for valid simulations only
### better than imported 'prox', which doesn't check time
p = (AllSums$prox==1) & (AllSums$t >= 1e7) & (AllSums$aBi > 0)
ProxSums = AllSums[p,]
attach(ProxSums)
#-------------------------------------------------------------------------
#skipped stuff
#-------------------------------------------------------------------------

###############################################################################
# periapse of original (aka binary-only) orbit
#pB2=((AllSums$aB+AllSums$da)*(1-(AllSums$eB+AllSums$de)))[p]
# periapse of final orbit
#pB3=AllSums$pB[p]

# a of 2- and 3-star system disks
#aB2 = (aB+da)[p]
#aB3 = aB[p]

# rtr of 2- and 3-star system disks
#rtr2 = r2[p]
#rtr3 = r3[p]

# compare median, mean, stdev, and rel stdev of rtr vs. a and rtr vs. p
r.a = cbind( c(median(r2/aBi),median(r3/aBf)), 
			 c(  mean(r2/aBi), mean(r3/aBf)), 
			 c(    sd(r2/aBi),   sd(r3/aBf)) )
	rownames(r.a) = c('bin','tri')
	colnames(r.a) = c('med','mean','sd')
relsd.a = r.a[,'sd']/r.a[,'mean']
r.a=cbind(r.a,relsd.a)
r.p = cbind( c(median(r2/pBi),median(r3/pBf)), 
			 c(  mean(r2/pBi),  mean(r3/pBf)), 
			 c(    sd(r2/pBi),    sd(r3/pBf)) )
	rownames(r.p) = c('bin','tri')
	colnames(r.p) = c('med','mean','sd')
relsd.p = r.p[,'sd']/r.p[,'mean']
r.p = cbind(r.p,relsd.p)


###############################################################################
### Compare rtr vs. p for B-2 vs B-3 -- normalized to aBin

f2=lm(r2~pBi)
f3=lm(r3~pBf)

for (iter in 1:2)	{
if (iter==1)	{
pdf('../Paper/Inserts/rTr-vs-peri.pdf',height=4.5,width=4.5)
	} else {
png('../Paper/Inserts/rTr-vs-peri.png',height=4.5,width=4.5, units="in",res=150)
	}
plot(pBi,r2, pch=20,
	xlab='Pericenter (AU)',
	ylab=expression('Truncation radius r'[tr]*' (AU)'),
	xlim=c(min(c( pBi, pBf)),max(c( pBi, pBf))),
	ylim=c(min(c(r2,r3)),
           max(c(r2,r3))) )
abline(f2)
points(pBf,r3, pch=20,col='red')
abline(f3,col='red')
legend('topleft',legend=c('Binary','Triple'),
	lty=1,pch=20,col=c('black','red'))
dev.off()
	} # iter

