### Indices for the outcomes of simulations
# ejection indices (single-B, single-C, single-either, double)
singB = ( !is.na(destB) & destB=='ejected' &  is.na(destC) )
singC = (  is.na(destB) & destC=='ejected' & !is.na(destC) )
sing  = (singB | singC)
doub  = !is.na(destB) & !is.na(destC)
# too huge (falsely counted as ejection in Mercury)
huge = ( ( (!is.na(destB) & eBf < 1.)   | 
		   (!is.na(destC) & eCf < 1.) ) & !is.na(EBf) & !is.na(ECf))
	huge[is.na(huge)]=FALSE
# collision, i.e. not 'ejected' and not NA/stable
coll  = (( !is.na(destB) & destB!='ejected' ) | 
		 ( !is.na(destC) & destC!='ejected' ) )
# error in reading => na values
brkn = is.na(aBf) | is.na(aCf) | is.na(eBf) | is.na(eCf)
# survival index
surv  = (is.na(destB) & is.na(destC) & eCf <1. & !brkn)
# growth index (did C move outward?)
grow  = (surv  & aCf>aC)
# proxima-like indice (10-20 kAU, 20+ kAU)
prox1 = (surv  & apCf>=r.prx & apCf<=r.big)
prox2 = (surv  & apCf>r.big)
prox  = (prox1 | prox2)

### Collect into one data frame
indices=data.frame(	surv,grow,prox1,prox2,prox,huge,
					singB,singC,sing,doub, coll, brkn)

### Data frame of plot point parameters
#cat=c('brkn','doub','singB','singC','coll', 'surv','grow','prox1','prox2','huge')
#colorRampPalette(brewer.pal(6,"Blues"))[6:2]	#(length(table.data))
#colorRampPalette(brewer.pal(6,"Reds" ))[2:6]

fate=rep( NA, length(surv))
	fate[surv]  ='surv'
	fate[grow]  ='grow'
	fate[prox1] ='prox1'
	fate[prox2] ='prox2'
	fate[huge]  ='huge'
	fate[doub]  ='doub'
	fate[coll]  ='coll'
	fate[singB] ='singB'
	fate[singC] ='singC'
fate=as.factor(fate)
pt=data.frame(	fate,
				cols=rep( NA, length(surv)), 
				pchs=rep( NA, length(surv)) )
	pt$pchs[surv|grow|prox1|prox2|huge]=19
	pt$pchs[doub|coll|singB|singC]     =4
	pt$pchs[brkn]=8	# alternately use separate print command to add pch 4
	pt$cols[surv]      ='blue'
	pt$cols[grow]      ='green'
	pt$cols[prox1|coll]='black'
	pt$cols[prox2]     ='grey'
	pt$cols[huge]      ='purple'
	pt$cols[doub]      ='red'
	pt$cols[singB]     ='orange'
	pt$cols[singC]     =colors()[144]

###############################################################################
### Write totals

### Correct for missing sims
old   = logtB>8 | logtC>8
early = Case=='Early'

# fate of sims that got past the 1e8 step:
rates = rbind( colMeans(indices[old & !early,]),
			   colMeans(indices[old &  early,]) )
	rownames(rates) = c('Late','Early')
corrections = rates*nmissing

### Get percent and number of each type
nums = rbind( colSums(indices[!early,]),
			  colSums(indices[ early,]) ) + corrections
	rownames(nums) = c('Late','Early')
pers = nums*100/(nsims+nmissing)
#s.long=data.frame(t(          (colSums(indices)+corrections)         ))
#m=data.frame(t( signif( (colMeans(indices)+corrections)*100, 3) ))
#m.long=s.long*100/(nsims+nmissing)
s=signif(t(nums), 3)
m=signif(t(pers), 3)

### print output
#cat(paste(nsims+nmissing,'recent simulations\n'))
#cat(paste("sims with no ejection  = ",m$surv ,'% (',s$surv ,')\n',sep=''))
#cat(paste("sims where C moves out = ",m$grow ,'% (',s$grow ,')\n',sep=''))
#cat(paste("Proxima-like sims      = ",m$prox ,'% (',s$prox ,')\n',sep=''))
#cat('\n')
#cat(paste("only B ejected         = ",m$singB,'% (',s$singB,')\n',sep=''))
#cat(paste("only C ejected         = ",m$singC,'% (',s$singC,')\n',sep=''))
#cat(paste("double ejection        = ",m$doub ,'% (',s$doub ,')\n',sep=''))
#cat(paste("Large orbit, false ejc = ",m$huge ,'% (',s$huge ,')\n',sep=''))
#cat('\n')
#cat(paste("Collision/other?       = ",m$coll ,'% (',s$coll ,')\n',sep=''))
#cat(paste("broken sims (NAs, etc) = ",m$brkn ,'% (',s$brkn ,')\n',sep=''))

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
br.e   =          1.*(0:n1)/n1
br.de  =  -1. +   2.*(0:n1)/n1
br.i   =         180*(0:n1)/n1
br.di  = -180.+2*180*(0:n1)/n1

### Plot cutoff line for proxima-like orbit (apocenter = 10,000 AU)
#d = 100.	# density of points in the line
a1.prx = round(r.prx/2):round(r.prx)
a2.prx = round(r.prx  ):round(r.prx*2)

e1.prx =   (r.prx/a1.prx)-1.
e2.prx = (2*r.prx/a2.prx)-1.

### Polygon vertices for shaded region showing most 'proxima-like' sims
vertices.e = c(   1.0,   1.0,  -0.0,  -0.0)
vertices.ap= c( 10000, 20000, 20000, 10000)
vertices.a = vertices.ap/(1+vertices.e)

