### Make plots and such for the suite of simulations indicated.
### run from command line:
#   R CMD BATCH -1 plot.R

### Get run version from input variable
args <- commandArgs(trailingOnly = F)
if (length(args) >2) options(echo = FALSE)
version <- sub("-","",args[length(args)])

print(args)

### If version not specified, try using version = 1
if (length(version) == 0 | 
	version=="-no-readline" | 
	version=="/usr/lib64/R/bin/exec/R") {
		print('no args')
		version = 1	}

### Packages
library(xtable)

### Settings
r.prx 	= 10000.		# min apoapse to be called 'proxima-like'
r.big 	= 20000.		# max apoapse to be called 'proxima-like'
DoCalcRunTime = TRUE	# Whether to run the calc time script

### Files to read data from
### Regular (current) version:
if (version==1)	{
	prefix='Plots/'			# Prefix for plot files
	sumfiles = c('SumAll.out')
	} else if (version==2)	{
### Older versions of sims (eps instead of E, no dE, unequal masses)
	prefix='Saved/CurrentMasses/'			# Prefix for plot files
	sumfiles=c(	'Saved/CurrentMasses/SumAll081414.out',
				'Saved/CurrentMasses/SumAll080614.out',
				'Saved/CurrentMasses/SumAll072314.out' )
	} else print('invalid version')

###############################################################################
### Calculate average runtime per sim
if (DoCalcRunTime) source('CalcRunTime.R')

###############################################################################
### Read in Summary.out, normal
source('ReadSumFiles.R')

###############################################################################
### Indices for the outcomes of simulations, with 2nd set including old sims
# error in reading => na values
brkn = is.na(aBf) | is.na(aCf) | is.na(eBf) | is.na(eCf)
# survival index
surv  = (is.na(destB) & is.na(destC) & eCf <1. & !brkn)
# growth index (did C move outward?)
grow  = (surv  & aCf>aC)
# proxima-like index
prox  = (surv  & aCf*(1+eCf)>=r.prx & aCf*(1+eCf)<=r.big)
# double-ejection index
doub  = !is.na(destB) & !is.na(destC)

# too huge (falsely counted as ejection in Mercury)
huge = ((!is.na(destB) & eBf < 1.) | (!is.na(destC) & eCf < 1.))
	huge[is.na(huge)]=FALSE

###############################################################################
### Write totals
i=signif(c(
	mean(brkn),mean(surv),mean(grow),mean(prox),mean(doub),mean(huge)
			)*100,3)
cat(paste(dim(SumAll)[1],'recent simulations\n'))
cat(paste("broken sims (NAs, etc) = ",i[1],'% (',sum(brkn),')\n',sep=''))
cat(paste("sims with no ejection  = ",i[2],'% (',sum(surv),')\n',sep=''))
cat(paste("sims where C moves out = ",i[3],'% (',sum(grow),')\n',sep=''))
cat(paste("Proxima-like sims      = ",i[4],'% (',sum(prox),')\n',sep=''))
cat(paste("double ejection        = ",i[5],'% (',sum(doub),')\n',sep=''))
cat(paste("Large orbit, false ejc = ",i[6],'% (',sum(huge),')\n',sep=''))

summaryrow=data.frame( nsims=dim(SumAll)[1], 
						surv=i[2], prox=i[4], doub=i[5], brkn=i[1] )
if (version==2)	rownames(summaryrow)=c('Current masses') else {
				rownames(summaryrow)=c('Equal masses')	}

print(xtable(summaryrow), file=paste(prefix,'summaryrow.tex',sep=''),
	only.contents    = TRUE,
	include.colnames = FALSE,
	hline.after      = NULL  )

###############################################################################
# Define cut and max times
tcut=1e5
tmax=max(cbind(tB,tC))
tmin=min(cbind(tB,tC))
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

###############################################################################
# change in B?
bgrowth=abs(aB-aBf)/aB
bgrowth[aBf<0]=NA

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
aBchange = rbind(aB[surv],aBf[surv],(aBf[surv]-aB[surv])/aB[surv]*100)
aCchange = rbind(aC[surv],aCf[surv],(aCf[surv]-aC[surv])/aC[surv]*100)

eBchange = rbind(eB[surv],eBf[surv],(eBf[surv]-eB[surv])/eB[surv]*100)
eCchange = rbind(eC[surv],eCf[surv],(eCf[surv]-eC[surv])/eC[surv]*100)

changes=rbind(round(aBchange[3,]),round(aCchange[3,]), round(eBchange[3,]),
	round(eCchange[3,]) ) #, round(iBchange[3,]),round(iCchange[3,]))
row.names(changes)=c('aB','aC','eB','eC') #,'iB*','iC*')
#print(changes)

###############################################################################
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

### Plot cutoff line for proxima-like orbit (apocenter = 10,000 AU)
d = 100.	# density of points in the line
a.prx = round(r.prx/2):round(r.prx)

e.prx = (r.prx/a.prx)-1.

###############################################################################
source('MakePlots.R')

###############################################################################
### Print proxima-like systems to a file (aCf>acutoff)
acutoff=5000	# AU
sink(paste(prefix,'proxlike.txt',sep=''))
options(width=300)
print(SumAll[prox,],row.names=FALSE)
options(width=80)
sink()

###############################################################################
### Print latex table of binary parameters in proxima-like systems

Bparams = data.frame( rep(' ',length(aB)),aB,eB,iB, aBf, eBf, iBf)[prox,]
	colnames(Bparams)[1] = ' '
print(xtable(Bparams), file=paste(prefix,'Bparams.tex',sep=''),
	only.contents    = TRUE,
	include.rownames = FALSE,
	include.colnames = FALSE,
	hline.after      = NULL  )

Cparams = data.frame( rep(' ',length(aB)), aC,eC,iC, aCf, eCf, iCf)[prox,]
	colnames(Cparams)[1] = ' '
print(xtable(Cparams), file=paste(prefix,'Cparams.tex',sep=''),
	only.contents    = TRUE,
	include.rownames = FALSE,
	include.colnames = FALSE,
	hline.after      = NULL  )

###############################################################################
detach(SumAll)

if (length(args) > 2) q('no')


