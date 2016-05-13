
#dir='C22/'
aeidir=paste(dir,'Aei/',sep='')
cent='AlCenB'

#
mSun   = 1.9891e33	# g
mEarth = 5.972e27	# g

#-----------------------------------------------------------------------------#
### read in number of objects and their names
if (dir=='C01/') {
#	n=3
#	files=c('AlCenA.aei','P0018.aei','P0053.aei') 
	n=(length(readLines(paste(dir,'/In/big.in',sep='')))-6)/4
	files=list.files(aeidir)
} else if (dir=='C06/') {
	n=3
	files=c('AlCenA.aei','P0003.aei','P0911.aei')
} else {
	n=(length(readLines(paste(dir,'/In/big.in',sep='')))-6)/4
	files=list.files(aeidir)
	}
print(paste(dir,': Nobj =',n,', nfiles =',length(files)))
names=sapply(strsplit(files, split='.', fixed=TRUE), function(x) (x[1]))

#-----------------------------------------------------------------------------#
### How many digits in obj names?
charlen=nchar(names[!names %in% c('AlCenA','AlCenB')][1])
Pcode=paste('P[[:digit:]]{',charlen-1,'}',sep='')

#-----------------------------------------------------------------------------#
### Find longest-surviving object
size=vector()
for (i in 1:length(names)) size[i]=file.info(paste(aeidir,files[i],sep=''))[1,1]
maxTind=(1:length(size))[size==max(size)][1]

### Find number of dimensions and length of sim
testdata = read.table(paste(aeidir,files[maxTind],sep=''),skip=4, header=F ) 
ndims = dim(testdata)[2]
if (ndims==12) {
	form='long'
	c.all=c('t','a','e','i','mass','dens','x','y','z', 'vx','vy','vz') 
	c.want=c('t','a','e','i','x','y','z','vx','vy','vz','mass')
} else { if (ndims==8) {
	form='short'
	c.all=c('t','a','e','i','peri','node','M','mass') 
	c.want=c('t','a','e','i','mass')
} }

colnames(testdata)=c.all
maxTlength = length(testdata$t)
minT       = testdata$t[2]
maxT       = max(testdata$t)

#-----------------------------------------------------------------------------#
### create array for data
# array of all possible numbers of stars that might appear in info files, 
# to replace them with NAs
na=c('*','**','***','****','*****','******','*******',
	'********','*********','**********','***********')

if (form=='short') {
aei = array(data = 0., dim = c( length(files), maxTlength, length(c.want)+2), 
	dimnames = list(names,NULL,c(c.want,'p','ap')) )
 } else if (form=='long') {
aei = array(data = 0., dim = c( length(files), maxTlength, length(c.want)+4), 
	dimnames = list(names,NULL,c(c.want,'p','ap','r','v')) )
 }

for (i in 1:length(names)) {
	filename=paste(aeidir,files[i],sep='')
	Tlength = length(readLines(filename) )-4
	objdata = read.table(filename, header=F,skip=4,	
						 col.names=c.all, na.strings=na)[,c.want]
	# It's stupid that I have to iterate over this, but otherwise it
	# turns my matrix into a list of n*maxTlength*ndim single elements...???
	for (j in 1:length(c.want))	{
		for (k in 1:Tlength)	aei[i, k, j] = objdata[ k,j]	}
	aei[i,,'p']  = aei[i,,'a']*(1-aei[i,,'e'])
	aei[i,,'ap'] = aei[i,,'a']*(1+aei[i,,'e'])
	if (form=='long')	{
		aei[i,,'r'] = (aei[i,, 'x']^2+aei[i,, 'y']^2+aei[i,, 'z']^2)^.5
		aei[i,,'v'] = (aei[i,,'vx']^2+aei[i,,'vy']^2+aei[i,,'vz']^2)^.5
		}	
	}

### Objects that survive to the end of the sim
Surv  = aei[,maxTlength,'t'] > 0
print(paste('Objects remaining:',sum(Surv)))
if (sum(Surv) < 50) {
	sink(paste(dir,'PlList.txt',sep=''))
	print(names[Surv])
	print(cbind(aei[Surv,maxTlength,c('a','e','i')],
	        m=aei[Surv,maxTlength,c('mass')]*mSun/mEarth) )
	sink()
}

### index for object that is the secondary star
isstar = names %in% c('AlCenA','AlCenB')

#-----------------------------------------------------------------------------#
### Determine fate of each object and assign plot color based on it
# skip lines starting with these words
ignorekeys = c('Fractional','Integration','Continuing','WARNING:','Modify',
	'Beginning','----------------------','-------------------','Algorithm:',
	'Output','Initial','Accuracy','Central',
	'J_2:','J_4:','J_6:','Ejection','Radius','Includes','Number','Integrating',
	'ERROR:', 'Check')
# get full info.out file
info = read.table(paste(dir,'Out/info.out',sep=''), sep='\n',stringsAsFactors=F)
# contents of info.out will be parsed into this data frame
interactions = data.frame(  obj  = character(),
							fate = character(),
							time = double(), stringsAsFactors=FALSE)
count = 0
for (i in 1:dim(info)[1])	{
	words = strsplit(info[i,],split='\\s+')[[1]][-1]
	if (!(words[1] %in% ignorekeys)) {
		count = count+1
		if        (length(words)==8) {interactions[count,]=c(words[5],words[1],words[7])
		} else if (length(words)==5) {interactions[count,]=c(words[1],'ejectd',words[4])
		} else if (length(words)==9) {interactions[count,]=c(words[1],    cent,words[8])
		} else {
			print(info[i,])
			stop('ERROR: unrecognized line in info.out file:')
		} # 'if # words' clauses
	} # 'if word not ignored' clause
} # i, lines of info
interactions$time = as.numeric(interactions$time)

### Find numer of objects between each colliding pair
dn = array(dim=dim(interactions)[1])
for (i in 1:length(dn))	{
	if ( length(grep(Pcode,interactions[i,'fate']))>0) {
		this = as.numeric(sub('P','',interactions[i, 'obj']))
		that = as.numeric(sub('P','',interactions[i,'fate']))
		dn[i] = this-that
	}
}
interactions$dn = dn

### Match each object to its fate from info file
fate = rep('remains',length(names))
for (i in 1:length(names))	{
	if (names[i] %in% interactions$obj)	{
		fate[i] = interactions$fate[interactions$obj == names[i]]
	}
}
#fate=as.factor(fate)

fcol = rep('black',length(names))
for (i in 1:length(fcol))	{
	if (fate[i]=='ejectd') fcol[i] = 'gray'
	if (fate[i]=='AlCenA') fcol[i] = 'orange'
	if (fate[i]=='AlCenB') fcol[i] = 'red'
	if (fate[i]=='remains') fcol[i] = 'blue'
	if ( length(grep(Pcode,fate[i]))>0) fcol[i] = 'lightblue'
	}

#-----------------------------------------------------------------------------#
### Ultimate fate (where did the mass end up?)
dest=c('AlCenA','AlCenB','ejectd','remains',names[Surv])
ult = rep('',length(fate))
ult[fate %in% dest] = fate[fate %in% dest]
ult[Surv] = names[Surv]

	while ( sum(ult=='')>0 )	{
		print(sum(ult==''))
		finalized = names[ult %in% dest]
		for (i in 1:length(finalized))	{
			ult[fate==finalized[i]] = ult[names==finalized[i]]	}
	}
allfates=levels(as.factor(ult))
starfate = ult[isstar]
if (length(starfate)==0) {
	print('Weirdness: non-binary simulation') 
 } else if (starfate != 'remains') {
	print(paste('Weirdness! Star fate =',starfate))
 }

#ult=ult[!names %in% c('AlCenA','AlCenB')]

### Calculate percentage of initial mass that meets each fate
nobj.fates=summary(as.factor(ult))
nobj.fracs=nobj.fates/sum(nobj.fates)
mass.fates=nobj.fates
	for (i in 1:length(nobj.fates)) {
		mass.fates[i]=sum(aei[
		(ult==attributes(nobj.fates)[[1]][i]) & !isstar,
		1,'mass']) }
mass.fracs = mass.fates/sum(mass.fates)

options(scipen=4)
sink(paste(dir,'FateFracs.txt',sep=''))
print('Fraction of mass meeting each fate:')
print(mass.fracs,digits=2)
sink()
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#### Function to plot the specified column vs. time
#vsT = function(ycol,lg='x'){

#	if ('AlCenA' %in% dimnames(aei)[[1]])	{
#		ind=which(dimnames(aei)[[1]]=='AlCenA')	
#	} else if ('AlCenB' %in% dimnames(aei)[[1]])	{
#		ind=which(dimnames(aei)[[1]]=='AlCenB')
#	} else {
#		ind=length(names)
#	}

#	xlims = c(minT,maxT)
#	if (ycol == 'e') {
#		ylims = c(0,1) 
#	} else if (ycol == 'i') {
#		ylims = c(0,180)
#	} else if (ycol == 'v') {
#		ylims = c(0,max(aei[2,,ycol]))
#	} else {
#		ylims = c(0, max(aei[ind,-1,ycol],na.rm=T))
#		}

#	lines= aei[ind,,'t'] > 0.

##	plot(aei[ind,lines,'t'], aei[ind,lines,ycol], 
##	log=lg, xlim=xlims, ylim=ylims, xlab='',ylab=ycol,xaxt='n',
##	pch=20,col=fcol[ind])
#	
##	counts=sort(summary(as.factor(fcol)),decreasing=T)
##	for (j in 1:length(counts)){
##		for (i in 1:length(files))	{
##			if (fcol[i] == attributes(counts)[[1]][j])	{
##				t = aei[i,-1,'t']
##				y = aei[i,-1,ycol]
##				e = aei[i,-1,'e']
##				points(t[e<=1.0], y[e<=1.0],pch=20,cex=.2,col=fcol[i])
##			}
##		}
##	}
#}
##-----------------------------------------------------------------------------#
#### Plot semimajor axis, eccentricity, pericenter, and apocenter vs. time
##pdf(paste(dir,'Tchanges.pdf',sep=''),width=10,height=10)
##png(paste(dir,'Tchanges.png',sep=''), res = 200, width = 10, height = 10, units = "in")
##par(mfrow=c(5,1), oma=c(5,0,1,0), mar=c(0,4,0,1))
##vsT('a')
##vsT('e')
##vsT('i')
##vsT('p')
##vsT('ap')
##axis(1)
##dev.off()
##-----------------------------------------------------------------------------#

#### Plot a-e
##pdf(paste(dir,'a-e.pdf',sep=''),width=10,height=10)
##plot(aei[1,,'a'], aei[1,,'e'], 
##	xlim=c( 0., 10. ), 
##	ylim=c( 0., 1.), 
##	pch=20,col=fcol[1])

##for (i in 2:length(files))	{
##	if (fcol[i] == 'gray')	{
##		a = aei[i,-1,'a']
##		e = aei[i,-1,'e']
##		points( a[e<=1.], e[e<=1.],pch=20,cex=.2, col=fcol[i])
##	}}
##for (i in 2:length(files))	{
##	if (fcol[i] != 'gray')	{
##		a = aei[i,-1,'a']
##		e = aei[i,-1,'e']
##		points( a[e<=1.], e[e<=1.],pch=20,cex=.2, col=fcol[i])
##	}}
##dev.off()

##-----------------------------------------------------------------------------#
##Function to make barplot of a parameter showing fates of objs
##library(plotrix)
##fatebars = function(xcol)	{
##	xmin = min(aei[ !isstar, 1, xcol ], na.rm=T)
##	xmax = max(aei[ !isstar, 1, xcol ], na.rm=T)
##	if (xcol!='mass') {
##		nbr  = 30
##		brk  = xmin + (0:nbr)*((xmax-xmin)/nbr)
##	} else {
##		lvls = as.numeric(levels(as.factor(aei[ !isstar, 1, xcol ])))
##		nbr = length(lvls)
##		brk = c(lvls[1]*.99,lvls*1.01)
##	}

##	h=array( dim=c( length(allfates), nbr ),dimnames=list(allfates,NULL) )
##	for (j in 1:length(allfates))	{
###		h[j,]=hist( aei[ ((ult==allfates[j]) & !isstar), 1, xcol ],
###			br=brk, plot=F)$counts
##		h[j,]=weighted.hist( x=aei[ ((ult==allfates[j]) & !isstar), 1, xcol ],
##			w=aei[ ((ult==allfates[j]) & !isstar), 1, 'mass' ],	
##			br=brk, plot=F)$counts
##		}
##	barplot(h, space=0, xlab=xcol, legend.text=T, 
##		col=c('black',rainbow(length(allfates)-1)), 
##		args.legend=list(x='top',ncol=2) )
##	axvals = axisTicks( c(xmin,xmax), log=F)
##	axlocs = (nbr+1)*(axvals-xmin)/(xmax-xmin)
##	axis(side=1,at=axlocs,labels=axvals)
###	return(h)
##	}
##fatebars('a')
##-----------------------------------------------------------------------------#
#### Where mass from initial disk ends up

#plotprms = c('a','mass')

#n=ceiling(sqrt(length(plotprms)))
#if (length(plotprms) <= n*(n-1)) m=n-1 else m=n

##pdf(paste(dir,'massfates.pdf',sep=''),height=m*6,width=n*6)
##par(mfrow=c(m,n), mar=c(5.1,2.1,1.1,1.1))

##for (i in 1:length(plotprms))	{
##	print(plotprms[i])
##	fatebars(plotprms[i])
##	}
##	
##dev.off()
###-----------------------------------------------------------------------------#
#### Distance from center and velocity vs. time
##if (form=='long')	{
##	png(paste(dir,'CartData.png',sep=''), res = 200, width = 10, height = 5, units = "in")
##	par(mfrow=c(2,1), oma=c(5,0,1,0), mar=c(0,4,0,1))
##	vsT('r',lg='')
##	vsT('v',lg='')
##	axis(1)
##	dev.off()
##	}
##-----------------------------------------------------------------------------#
#### Get track of star over time
#if ('AlCenA' %in% names)	{

#A=aei['AlCenA',,]
#Ais= ((A[,'a'] != 0.) | is.na(A[,'a']))

##pdf(paste(dir,'AlCenApos.pdf',sep=''), width=10, height=10)
##par(mfrow=c(2,2))

##plot(A[Ais,'x'],A[Ais,'y'],type='l')
##points(A[Ais,'x'],A[Ais,'y'],pch=20,col='blue')
##points(0.,0., pch=19, col='orange')

##plot(A[Ais,'x'],A[Ais,'z'],type='l')
##points(A[Ais,'x'],A[Ais,'z'],pch=20,col='blue')
##points(0.,0., pch=19, col='orange')

##plot(A[Ais,'y'],A[Ais,'z'],type='l')
##points(A[Ais,'y'],A[Ais,'z'],pch=20,col='blue')
##points(0.,0., pch=19, col='orange')

##dev.off()
##-----------------------------------------------------------------------------#

#r=A[Ais,'r']
#t=A[Ais,'t']*365.25
#v=A[Ais,'v']

##r=r[5:length(r)]
##t=t[5:length(t)]
##v=v[5:length(v)]


#dr=array()
#dt=array()
#dv=array()

#for (i in 1:(length(r)-1))	{
#	dr[i]=r[i+1]-r[i]
#	dt[i]=t[i+1]-t[i]
#	dv[i]=v[i+1]-v[i]
#	}

#drt=abs(dr/dt)
#dv=abs(dv)


#up=(abs(dt-4.)<.1)
#dn=(abs(dt-3.)<.1)


#if (dir %in% c('C01/','C06/','C26/','C28/')) {

#if (dir=='C01/') {
#xl=c(5.21e4,5.28e4)
#yl1=c(0,20)
#yl2=c(0,0.05)
#} else if (dir=='C06/') {
#xl=c(2.294e5,2.34e5)
#yl1=c(0,100)
#yl2=c(0,0.03)
#} else if (dir=='C26/') {
#xl=c(365.25e3,1200000)	# t scale
#yl1=c(0,50)			# r scale
#yl2=c(0,0.01)			# v scale
#} else if (dir=='C28/') {
#xl=c(.9e6,1.1e6)
#yl1=c(0,50)			# r scale
#yl2=c(0,0.01)
#} else if (dir=='C11/') {
#xl=c(2.5e4,3e4)
#yl1=c(0,75)			# r scale
#yl2=c(0,0.04)
#}

##pdf(paste(dir,'weirdness.pdf',sep=''),width=10,height=10)
##par(mfrow=c(4,1))
##plot(    t,   r, pch=20, xlim=xl, ylim=yl1)
##plot(    t,   v, pch=20, xlim=xl, ylim=yl2)
##plot(t[-1], drt, pch=20, xlim=xl, ylim=yl2)
##plot(t[-1], drt, pch=20, xlim=xl)
##dev.off()

#	}	# if C01 or 6
#}	#if star exists

