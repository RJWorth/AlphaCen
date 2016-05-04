
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
#dirs=c('Bin03','Bin04','Bin09','Bin10',
#       'Bin01','Bin02','Bin07','Bin08',
#       'Bin05','Bin06','Bin11','Bin12',
#	   'Bin15','Bin13','Bin16','Bin14',
#       'Bin17','Bin19','Bin18','Bin20',
#       'Bin21','Bin23','Bin22','Bin24')

### Physical properties of the above simulations
#	alpha = c(rep(1.5,12), rep(1.,12))
#	sig = rep( c(rep(3,4), rep(1,4), rep(.3,4)), 2 )
#	rtr = rep(c(2.77, 3.08),12)
#	rej = rep( 1e2, 24)

### Read in sim names and parameters
simparams = read.table('BinSimCatalogue.txt',header=T, as.is=T)
simparams = simparams[1:24,] # they don't all exist yet
dirs  = simparams$dir
alpha = simparams$alpha
sig   = simparams$sig
rtr   = simparams$rtr

### Locations
supdir = 'Data/'
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
simparams[,'ej.frac'] = ej.frac
simparams[,'st.frac'] = st.frac
simparams[,'nplts']   = nplts
simparams[,'n3big']   = n3big

simparams$dir   = as.factor(simparams$dir)
simparams$sig   = as.factor(simparams$sig)
simparams$alpha = as.factor(simparams$alpha)
simparams$rtr   = as.factor(simparams$rtr)

plot(simparams$ej.frac, simparams$n3big, lwd=2,
    pch = as.numeric(simparams$alpha) +20, 
    col = as.numeric(simparams$rtr  ) , 
    bg  = as.numeric(simparams$sig  ) )


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
### Read in the functions that make the plots
source('PltFns.R')
### Plot planet systems
# one plot for each rtr value
for (j in levels(simparams$rtr)) {
	plottag = sub('.','_',j, fixed=T)
	pdf(paste(tag,'PltSystems-',plottag,'.pdf',sep=''), width=6.5, height=5.)
	nplots = length(which(simparams$rtr == j))
	par(mfrow=c((nplots+2)/2,2), oma=c(5,5,1,2))
	# iterate through each sigma level
	for (k in sort(levels(simparams$sig),decreasing=T)) {
		# alternate alpha levels (L/R on plot)
		al1 = which((simparams$rtr == j) & (simparams$sig == k) & (simparams$alpha==1))
		al2 = which((simparams$rtr == j) & (simparams$sig == k) & (simparams$alpha==1.5))
		stopifnot(length(al1)==length(al2))
		for (l in 1:length(al1)) {
			# plot left panel for this row
			print(simparams[al1[l],])
			writeax = F                       # don't plot x axis labels
			par(mar=c(0,1,0,0))               # no space between plots
			PlotSim(al1[l], 'L', writeax)
			# plot right panel
			print(simparams[al2[l],])
			par(mar=c(0,0,0,1))
			if ((k==levels(simparams$sig)[1]) & (l == length(al1))) writeax = T
			PlotSim(al2[l], 'R', writeax)
		} #l gives line in simparams to use (alternating alpha)
	} # k, sig
	PlotSS()
	PlotLegend()
	MarginText()
	dev.off()
} # j, rtr

###############################################################################

