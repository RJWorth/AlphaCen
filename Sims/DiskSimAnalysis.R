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

### Index to select for valid simulations only
### better than imported 'prox', which doesn't check time
p = (AllSums$prox==1) & (AllSums$t >= 1e7) & (AllSums$aB > 0)
#-------------------------------------------------------------------------
#skipped stuff
#-------------------------------------------------------------------------


###############################################################################
### Compare rtr vs. p for B-2 vs B-3 -- normalized to aBin
pB2=((AllSums$aB+AllSums$da)*(1-(AllSums$eB+AllSums$de)))[p]
pB3=AllSums$pB[p]

f2=lm(AllSums$r2[p]~pB2)
f3=lm(AllSums$r3[p]~pB3)

for (iter in 1:2)	{
if (iter==1)	{
pdf('../Paper/Inserts/rTr-vs-peri.pdf',height=4.5,width=4.5)
	} else {
png('../Paper/Inserts/rTr-vs-peri.png',height=4.5,width=4.5, units="in",res=150)
	}
plot(pB2,AllSums$r2[p], pch=20,
	xlab='Pericenter (AU)',
	ylab=expression('Truncation radius r'[tr]*' (AU)'),
	xlim=c(min(c( pB2, pB3)),max(c( pB2, pB3))),
	ylim=c(min(c(AllSums$r2[p],AllSums$r3[p])),
           max(c(AllSums$r2[p],AllSums$r3[p]))) )
abline(f2)
points(pB3,AllSums$r3[p], pch=20,col='red')
abline(f3,col='red')
legend('topleft',legend=c('Binary','Triple'),
	lty=1,pch=20,col=c('black','red'))
dev.off()
	} # iter

