
#dirs=c('C05','C07','C08', 'C13','C14','C15','C16','C17','C18','C19','C20', 'C21','C22','C23','C24','C25','C27')
#dirs=c('C35','C36','C37')

### 'no A' directories
#dirs=c('C22','C23','C35','C29','C30','C36','C24','C25','C37')
#	rtr = c(2.77, 2.77, 2.77, 2.77, 2.77, 2.77, 3.08, 3.08, 3.08)
#	rej = c(2.77, 2.77, 2.77, 3.08, 3.08, 3.08, 3.08, 3.08, 3.08)
#tag='../Heuristic/NoA-'
#### 'with A' directories w/o discontinuities
#dirs=c('C21', 'C31', 'C32', 'C27', 'C33')
#	rtr = c(2.77, 2.77, 2.77, 3.08, 3.08, 2.77)
#	rej = c( 1e2,  1e2,  1e2,  1e2,  1e2, 1e2)
### New simulations with Binary integrator
dirs=c('Bin01','Bin02','Bin03','Bin04','Bin05','Bin06', 
	   'Bin07','Bin08','Bin09','Bin10','Bin11','Bin12')
supdir = 'Data/'
	rtr = rep(c(2.77, 3.08),6)
	sig = rep(c( 1, 1, 3, 3,.3,.3),2)
	rej = rep( 1e2, 12)

tag='../Paper/Inserts/'

elem=list()

###
mSun   = 1.9891e33	# g
mEarth = 5.972e27	# g

###
for (i in 1:length(dirs)){
	fname = paste('Data/',dirs[i],'/Out/element.out',sep='')
	elem[[i]] = read.table(fname,skip=3)
	elem[[i]][,'mass'] = elem[[i]][,'mass']*mSun/mEarth
	
	print(paste( dirs[i], ':',dim(elem[[i]])[1],dim(elem[[i]])[2]  ))
}
notstar = rownames(elem[[1]]) != 'AlCenA'
nplts = array(dim(elem[[1]][notstar,])[1])
Dir = rep(toString(1),nplts[1])
rTr = rep(     rtr[1], nplts[1])
sg  = rep(     sig[1], nplts[1])
theseplts = cbind(Dir, rTr, sg, elem[[1]][notstar,])
plts = theseplts
for (i in 2:length(dirs))	{
	notstar = rownames(elem[[i]]) != 'AlCenA'
	nplts[i] = dim(elem[[i]][notstar,])[1]
	Dir = rep(toString(i), nplts[i])
	rTr = rep(     rtr[i], nplts[i])
	sg  = rep(     sig[i], nplts[i])
	theseplts = cbind(Dir, rTr, sg, elem[[i]][notstar,])
	plts=rbind( plts, theseplts )
	}

massfracs = list()
for (ind in 1:length(dirs))	{
#	dir=paste(dirs[ind],'/',sep='')
	dir=paste(supdir,dirs[ind],'/',sep='')
	source('ReadDisk.R')
	print(ind)
	print(mass.fracs)
	massfracs[[ind]] = mass.fracs
	print(massfracs[[ind]])
	}

ej.frac = array()
for (i in 1:length(dirs)) if ('ejectd' %in% names(massfracs[[i]])) {
	ej.frac[i] = massfracs[[i]]['ejectd'] } else ej.frac[i] = 0.
print(ej.frac)


fates=c('ejectd','AlCenA','AlCenB',sapply(1:max(nplts), function(x) paste('P',x,sep='')))
Fates = array(dim=c(length(dirs),length(fates)), dimnames=list( (dirs),(fates) ) )

for (ind in 1:length(dirs))	{
	sortfracs = sort(massfracs[[ind]], decreasing=T)
	pltind = grep('P[[:digit:]]{1}',names(sortfracs))

	for (ind2 in 1:3) 				Fates[ind,ind2] = sortfracs[fates[ind2]]
	for (ind2 in 1:length(pltind)) {
		Fates[ind,ind2+3] = sortfracs[pltind[ind2]]
		}
	}

sink('PltSystemSummaries.txt')
options(width=200)
print(Fates)
options(width=80)
sink()

###############################################################################

### Text for which rtr and sigma values were used in each sim, to add to plots
SimLabs = vector("expression",length(dirs))
for (i in 1:length(dirs)) { 
#	quote1 = bquote( r[tr] == .(rtr[i])*',' ~ c[sigma] == .(sig[i]) )
	SimLabs[i] = substitute(expression( r[tr] == r. ~ 'AU,' ~ c[sigma] == s.),
					list(r. = rtr[i], s. = sig[i]) )[2]
}

### Index for which objects are stable planets
plind=list()
for (i in 1:length(dirs))	{
plind[[i]]=(grepl(paste('P[[:digit:]]{1,}',sep=''), rownames(elem[[i]]) ) & 
			elem[[i]][,'a']>0.)
}

### max and min semimajor axes of planets
alims = c(1,1)
for (i in 1:length(dirs))	{
	if (min(elem[[i]][plind[[i]],'a']) < alims[1]) alims[1]=min(elem[[i]][plind[[i]],'a'])
	if (max(elem[[i]][plind[[i]],'a']) > alims[2]) alims[2]=max(elem[[i]][plind[[i]],'a'])
	}
elims = c(-0.05,ceiling(max(plts$e)*10)/10+0.05)

#alims[2] = max(alims[2],rtr)

### Solar System planet parameters
SS = array(dim=c(4,6), 
	dimnames=list(c('Mercury','Venus','Earth','Mars'),colnames(elem[[1]])))
SS['Mercury',] = c(  0.387098,   0.205630,   7.005, 0.055, NA, NA)
SS['Venus',]   = c(  0.723332, 0.00677323, 3.39458, 0.815, NA, NA)
SS['Earth',]   = c(1.00000261, 0.01671123, 0.00005,    1., NA, NA)
SS['Mars',]    = c(  1.523679,     0.0935,   1.850, 0.107, NA, NA)

maxplt = array(dim=c(length(dirs)+1,2),
	dimnames=list(c(dirs,'SS'),c('mass','e') ))
for (i in 1:length(dirs)) maxplt[i,] = c(max(elem[[i]][plind[[i]],'mass']),
										 max(elem[[i]][plind[[i]],'e']))
	maxplt[length(dirs)+1,] = c(max(SS[,'mass']),max(SS[,'e']))

### Make a transparent gray color to use for HZ boxes
diskcol = rgb(190, 190, 190, alpha=100, maxColorValue=255)
HZcol   = rgb(176, 226, 255, alpha=100, maxColorValue=255)

### Reorder the dirs to plot like with like
ord = c(5,6,11,12,1,2,7,8,3,4,9,10)
###############################################################################
### Plot planet systems -- currently each system symbol size is normalized to its own max
### need to figure out how to make it consistent between systems
pdf(paste(tag,'PltSystems.pdf',sep=''), width=6, height=6.)
#png(paste(tag,'PltSystems.png',sep=''), width=6, height=6., units="in",res=900)
#setEPS(horizontal=F, onefile=F, paper='special')
#postscript(paste(tag,'PltSystems.eps',sep=''), width=6, height=6.)
### two columns
par(mfrow=c(ceiling((length(dirs)+1)/2),2), oma=c(5,5,5,2))
for (i in ord)	{
	writeax = F
	if (i %% 2 == 0) {
		par(mar=c(0,0,0,1))
		if (i == ord[length(ord)]) writeax = T
	} else {
		par(mar=c(0,1,0,0))
	}
	plot(0.,0., pch=NA, xlab='',ylab='',, xaxt='n', yaxt='n',
	xlim=c(0,alims[2]), ylim=elims)

	rect(xleft=0.693,xright=1.241, ybottom=-1, ytop=2, col=HZcol, lty=0)
	rect(xleft=0.35, xright=rtr[i],  ybottom=-0.01, ytop=0.01, col=diskcol, lty=0)

	points(elem[[i]][1:(nplts[i]-1),'a'], elem[[i]][1:(nplts[i]-1),'e'], pch=20, cex=3.*elem[[i]][1:(nplts[i]-1),'mass']**(1./2.))
	
	# X axis
	Axis(side=1, labels=writeax)
	# Y axis
	if (i %% 2 == 0) {
		axis(side=4, at = c(0, .1, .2, .3), labels=c('0','','0.2',''),cex.axis=.9)
	} else {
		Axis(side=2, at = c(0, .1, .2, .3), labels=c('0','','0.2',''),cex.axis=.9)
	}
	text(0,.9*elims[2], SimLabs[i], adj = c(0.,1))
	
	}
# SS planets
par(mar=c(0,1,0,0))
plot(0.,0., pch=NA, xlab='',ylab='',, xaxt='n', yaxt='n',
xlim=c(0,alims[2]), ylim=c(-0.05,ceiling(max(plts$e)*10)/10+0.05) )

rect(xleft=0.951,xright=1.676, ybottom=-1, ytop=2, col=HZcol, lty=0)

points(SS[,'a'], SS[,'e'], pch=20, cex=3.*SS[,'mass']**(1./2.))

Axis(side=1, labels=T)
Axis(side=2, at = c(0, .1, .2, .3), labels=c('0','','0.2',''),cex.axis=.9)
text(0.,.9*elims[2],'Solar System', adj = c(0.,1))

### For 2-column version, create empty plot to hold legend
plot(1, type="n", axes=F, xlab="", ylab="")
legend('bottom', legend=c('1 earth mass','2 earth masses','3 earth masses'), 
	pch=20, pt.cex=3.*c(1,2,3)**(1./2), bg='white')
### Outer margin modifications
mtext(side=1,line=3,outer=T,'a (AU)')
mtext(side=2,line=3,outer=T,'e')
mtext(side=3,line=1,outer=T,text='Planetary Systems Formed')

dev.off()
###############################################################################

