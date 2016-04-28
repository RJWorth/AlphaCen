
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
#dirs=c('Bin01','Bin02','Bin03','Bin04','Bin05','Bin06', 
#	   'Bin07','Bin08','Bin09','Bin10','Bin11','Bin12',
#	   'Bin13','Bin14','Bin15','Bin16','Bin17','Bin18', 
#	   'Bin19','Bin20','Bin21','Bin22','Bin23','Bin24')
dirs=c('Bin03','Bin04','Bin09','Bin10',
       'Bin01','Bin02','Bin07','Bin08',
       'Bin05','Bin06','Bin11','Bin12',
	   'Bin15','Bin13','Bin16','Bin14',
       'Bin17','Bin19','Bin18','Bin20',
       'Bin21','Bin23','Bin22','Bin24')
supdir = 'Data/'

### Physical properties of the above simulations
	alpha = c(rep(1.5,12), rep(1.,12))
	sig = rep( c(rep(3,4), rep(1,4), rep(.3,4)), 2 )
	rtr = rep(c(2.77, 3.08),12)
	rej = rep( 1e2, 24)

#tag='../Paper/Inserts/'
tag='../Paper/2ndDraft/'

elem=list()

### Constants
mSun    = 1.9891e33	# g
mEarth  = 5.972e27	# g
mPbig   = 3.22604695591e-07*mSun/mEarth # from mSun to mEarth units
mPsmall = 3.69396867427e-08*mSun/mEarth # from mSun to mEarth units
PlLim   = 3.*mPbig   # how many mMars an object must be to count as a planet

### Make data frame of each planet, with data about simulation it comes from
for (i in 1:length(dirs)){
	fname = paste('Data/',dirs[i],'/Out/element.out',sep='')
	elem[[i]] = read.table(fname,skip=3)
	elem[[i]][,'mass'] = elem[[i]][,'mass']*mSun/mEarth
	
	print(paste( dirs[i], ':',dim(elem[[i]])[1],dim(elem[[i]])[2]  ))
}
notstar = rownames(elem[[1]]) != 'AlCenA'
nplts = array(dim(elem[[1]][notstar,])[1])
Dir = rep( substr(dirs[1], nchar(dirs[1])-1, nchar(dirs[1])), nplts[1])
al  = rep(   alpha[1], nplts[1])
sg  = rep(     sig[1], nplts[1])
rTr = rep(     rtr[1], nplts[1])
theseplts = cbind(Dir, al, sg, rTr, elem[[1]][notstar,])
n1big = array(sum(theseplts$mass >= mPbig))
n2big = array(sum(theseplts$mass >= mPbig*2.))
n3big = array(sum(theseplts$mass >= mPbig*3.))
plts = theseplts
for (i in 2:length(dirs))	{
	notstar = rownames(elem[[i]]) != 'AlCenA'
	nplts[i] = dim(elem[[i]][notstar,])[1]
	Dir = rep( substr(dirs[i], nchar(dirs[i])-1, nchar(dirs[i])), nplts[i])
	al  = rep(   alpha[i], nplts[i])
	sg  = rep(     sig[i], nplts[i])
	rTr = rep(     rtr[i], nplts[i])
	theseplts = cbind(Dir, al, sg, rTr, elem[[i]][notstar,])
	n1big[i] = array(sum(theseplts$mass >= mPbig))
	n2big[i] = array(sum(theseplts$mass >= mPbig*2.))
	n3big[i] = array(sum(theseplts$mass >= mPbig*3.))
	plts=rbind( plts, theseplts )
	}

# up to here is fast, but this part takes a bit
### Get fraction of disk mass that goes to each planet, star, or ejection
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

### Ejected mass fraction in each sim
ej.frac = array()
for (i in 1:length(dirs)) if ('ejectd' %in% names(massfracs[[i]])) {
	ej.frac[i] = massfracs[[i]]['ejectd'] } else ej.frac[i] = 0.
### Mass fraction accreted onto stars in each sim
st.frac = array()
for (i in 1:length(dirs)) {
	if ('AlCenB' %in% names(massfracs[[i]])) {
	st.frac[i] = massfracs[[i]]['AlCenB'] } else st.frac[i] = 0.
	if ('AlCenA' %in% names(massfracs[[i]])) {
	st.frac[i] = massfracs[[i]]['AlCenA'] + st.frac[i] }
	}
### Mass fraction lost to either stars or ejection
lost.frac = ej.frac+st.frac

### dataframe of summary of each sim
sims = data.frame(cbind(as.numeric(substr(dirs, 4, 5)), alpha, sig, rtr, ej.frac, 
             nplts, n1big, n2big, n3big))
colnames(sims)[1]='Dir'
sims$Dir   = as.factor(sims$Dir)
sims$sig   = as.factor(sims$sig)
sims$alpha = as.factor(sims$alpha)
sims$rtr   = as.factor(sims$rtr)

plot(sims$ej.frac, sims$n1big, lwd=2,
    pch = as.numeric(sims$alpha) +20, 
    col = as.numeric(sims$rtr  ) , 
    bg  = as.numeric(sims$sig  ) )


### array of how much mass went to each destination
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
	SimLabs[i] = substitute(
        expression( r[tr] == r. ~ 'AU,' ~ alpha == a. ~ ',' ~ c[sigma] == s.),
        list(a. = alpha[i] , r. = rtr[i], s. = sig[i]))[2]
}

### Index for which objects are stable planets, and planetesimals
plind=list()
pmind=list()
for (i in 1:length(dirs))	{
plind[[i]]=(grepl(paste('P[[:digit:]]{1,}',sep=''), rownames(elem[[i]]) ) & 
			elem[[i]][,'a']>0. & 
            elem[[i]][,'mass'] > PlLim)
pmind[[i]]=(grepl(paste('P[[:digit:]]{1,}',sep=''), rownames(elem[[i]]) ) & 
			elem[[i]][,'a']>0. & 
            elem[[i]][,'mass'] < PlLim)
}

### max and min semimajor axes of planets
alims = c(1,1)
for (i in 1:length(dirs))	{
	if (min(elem[[i]][plind[[i]],'a']) < alims[1]) alims[1]=min(elem[[i]][plind[[i]],'a'])
	if (max(elem[[i]][plind[[i]],'a']) > alims[2]) alims[2]=max(elem[[i]][plind[[i]],'a'])
	}
alims = c(   0, alims[2]*1.1)
elims = c(-0.1,  ceiling(max(plts$e)*10)/10+0.1)

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
#ord = c(5,6,11,12,1,2,7,8,3,4,9,10)
ord=1:24

### Outer line thicknesses
thick = 13.
med   = 6.
### text position
textpos = c(0.1,.85*elims[2],alims[2]*.95)

###############################################################################
### Plot planet systems -- currently each system symbol size is normalized to its own max
### need to figure out how to make it consistent between systems
#png(paste(tag,'PltSystems.png',sep=''), width=6, height=6., units="in",res=900)
#setEPS(horizontal=F, onefile=F, paper='special')
#postscript(paste(tag,'PltSystems.eps',sep=''), width=6, height=6.)
pdf(paste(tag,'PltSystems.pdf',sep=''), width=7.5, height=9.)
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
	xlim=c(alims[1],alims[2]), ylim=elims, xaxs='i', yaxs='i')

	rect(xleft=0.693,xright=1.241, ybottom=-1, ytop=2, col=HZcol, lty=0)
	rect(xleft=0.35, xright=rtr[i],  ybottom=-0.01, ytop=0.01, col=diskcol, lty=0)
	### Plot planets (> some mass cutoff) black, smaller objs gray
	points(elem[[i]][plind[[i]],'a'], elem[[i]][plind[[i]],'e'], 
			pch=20, cex=3.*elem[[i]][plind[[i]],'mass']**(1./2.))
	points(elem[[i]][pmind[[i]],'a'], elem[[i]][pmind[[i]],'e'], 
			pch=20, cex=3.*elem[[i]][pmind[[i]],'mass']**(1./2.), col='gray')
	
	# X axis
	Axis(side=1, labels=writeax)
	# Y axis
	if (i %% 2 == 0) {
		axis(side=4, at = c(0, .1, .2, .3, .4, .5), las=2,
                                   labels=c('0','','0.2','','0.4',''),cex.axis=.9)
	} else {
		Axis(side=2, at = c(0, .1, .2, .3, .4, .5), las=2,
                                   labels=c('0','','0.2','','0.4',''),cex.axis=.9)
	}
	text(textpos[1],textpos[2], SimLabs[i], adj = c(0.,1))
	text(textpos[3],textpos[2], adj = c(1,1), 
		labels=paste(sprintf("%.1f", round(lost.frac[i]*100,1)),'% mass loss',sep=''))

	# add thicker lines to outer left and right sides
	if (i %% 2 == 1) {
		abline(v=alims[1], lw=thick)
		abline(v=alims[1], lw=1., col='red')
		abline(v=alims[2], lw=thick/2)
	} else if (i %% 2 == 0) {
		abline(v=alims[1], lw=thick/2)
		abline(v=alims[2], lw=thick)
		}
	# draw thicker lines between different sigma subsets
	if (i %in% c(1,2)) abline(h=elims[2],lw=thick)
	if (i %in% c(1,2, 5,6,  9,10, 13,14, 17,18, 21,22)) abline(h=elims[2],lw=med)
	if (i %in% c(23,24)) abline(h=elims[1],lw=thick)
	# draw an extra thick line between the alpha subsets
	if (i %in% c(11,12)) {
		abline(h=elims[1],lw=thick/2) }		
	if (i %in% c(13,14)) {
		abline(h=elims[2],lw=thick/2) }

	}
# SS planets
par(mar=c(0,1,0,0))
plot(0.,0., pch=NA, xlab='',ylab='',, xaxt='n', yaxt='n', xaxs='i', yaxs='i',
xlim=c(alims[1],alims[2]), ylim=c(-0.05,ceiling(max(plts$e)*10)/10+0.05) )

rect(xleft=0.951,xright=1.676, ybottom=-1, ytop=2, col=HZcol, lty=0)

points(SS[,'a'], SS[,'e'], pch=20, cex=3.*SS[,'mass']**(1./2.))

Axis(side=1, labels=T)
Axis(side=2, at = c(0, .1, .2, .3), labels=c('0','','0.2',''),cex.axis=.9,las=2)
text(textpos[1], textpos[2],'Solar System', adj = c(0.,1))

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

