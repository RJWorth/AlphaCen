
dir='C05/'
aeidir=paste(dir,'Aei/',sep='')
cent='B'

#-----------------------------------------------------------------------------#
### read in number of objects and their names
n=(length(readLines(paste(dir,'/In/big.in',sep='')))-6)/4 # = 1722
files=list.files(aeidir) # =1722
names=sapply(strsplit(files, split='.', fixed=TRUE), function(x) (x[1]))

#-----------------------------------------------------------------------------#
### Find longest-surviving object
size=vector()
for (i in 1:length(names)) size[i]=file.info(paste(aeidir,files[i],sep=''))[1,1]
maxTind=(1:length(size))[size==max(size)][1]

### Find number of dimensions and length of sim
testdata = read.table(paste(aeidir,files[maxTind],sep=''),skip=4, header=F ) 
ndims = dim(testdata)[2]
if (ndims==12) {
	c.all=c('t','a','e','i','mass','dens','x','y','z', 'vx','vy','vz') 
	c.want=c('t','a','e','i','x','y','z','vx','vy','vz','mass')
} else { if (ndims==8) {
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
	'********','*********','**********')

aei = array(data = 0., dim = c( length(files), maxTlength, length(c.want)+2), 
	dimnames = list(names,NULL,c(c.want,'p','ap')) )

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
	}

### 
isstar = names %in% c('AlCenA','AlCenB')

#-----------------------------------------------------------------------------#
### Determine fate of each object and assign plot color based on it
info = read.table(paste(dir,'Out/info.out',sep=''), sep='\n',stringsAsFactors=F)
interactions = data.frame(  obj  = character(),
							fate = character(),
							time = double(), stringsAsFactors=FALSE)

ignorekeys = c('Fractional','Integration','Continuing','WARNING:','Modify',
	'Beginning','----------------------','-------------------','Algorithm:',
	'Output','Initial','Accuracy','Central',
	'J_2:','J_4:','J_6:','Ejection','Radius','Includes','Number','Integrating')
count = 0
for (i in 1:dim(info)[1])	{
	words = strsplit(info[i,],split='\\s+')[[1]][-1]
	if (!(words[1] %in% ignorekeys)) {
		count = count+1
		if (length(words)==8) interactions[count,]=c(words[5],words[1],words[7])
		if (length(words)==9) interactions[count,]=c(words[1],'AlCenB',words[8])
		if (length(words)==5) interactions[count,]=c(words[1],'ejectd',words[4])
	}
}
interactions$time = as.numeric(interactions$time)

### Find numer of objects between each colliding pair
dn = array(dim=dim(interactions)[1])
for (i in 1:length(dn))	{
	if ( length(grep('P[[:digit:]]{4}',interactions[i,'fate']))>0) {
		this = as.numeric(sub('P','',interactions[i, 'obj']))
		that = as.numeric(sub('P','',interactions[i,'fate']))
		dn[i] = this-that
		print(paste(this,that,dn[i]))
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
	if ( length(grep('P[[:digit:]]{4}',fate[i]))>0) fcol[i] = 'lightblue'
	}

#-----------------------------------------------------------------------------#
### Ultimate fate (where did the mass end up?)
dest=c('AlCenA','AlCenB','ejectd','remains')
ult = fate
ult[grep('P[[:digit:]]{4}',fate)] = ''

	while ( sum(ult=='')>0 )	{
		print(sum(ult==''))
		finalized = names[ult %in% dest]
		for (i in 1:length(finalized))	{
			ult[fate==finalized[i]] = ult[names==finalized[i]]	}
	}
starfate = ult[isstar]
if (starfate != 'remains') print(paste('Weirdness! Star fate =',starfate))

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
print('Fraction of mass meeting each fate:')
print(mass.fracs)
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
### Function to plot the specified column vs. time
vsT = function(ycol){

	if ('AlCenA' %in% dimnames(aei)[[1]])	{
		ind=which(dimnames(aei)[[1]]=='AlCenA')	
	} else if ('AlCenB' %in% dimnames(aei)[[1]])	{
		ind=which(dimnames(aei)[[1]]=='AlCenB')
	} else {
		ind=length(names)
	}

	xlims = c(5e4,3e5)
	if (ycol == 'e') ylims = c(0,1) else if (ycol == 'i') ylims = c(0,180) else {
		ylims = c(0, max(aei[ind,-1,ycol],na.rm=T))	}

	plot(aei[ind,-1,'t'], aei[ind,-1,ycol], 
	log='x', xlim=xlims, ylim=ylims, xlab='',ylab=ycol,xaxt='n',
	pch=20,col=fcol[ind])
	
	counts=sort(summary(as.factor(fcol)),decreasing=T)
	for (j in 1:length(counts)){
		for (i in 1:length(files))	{
			if (fcol[i] == attributes(counts)[[1]][j])	{
				t = aei[i,-1,'t']
				y = aei[i,-1,ycol]
				e = aei[i,-1,'e']
				points(t[e<=1.0], y[e<=1.0],pch=20,cex=.2,col=fcol[i])
			}
		}
	}
}
#-----------------------------------------------------------------------------#
### Plot semimajor axis, eccentricity, pericenter, and apocenter vs. time
pdf(paste(dir,'Tchanges.pdf',sep=''),width=10,height=10)
par(mfrow=c(5,1), oma=c(5,0,1,0), mar=c(0,4,0,1))
vsT('a')
vsT('e')
vsT('i')
vsT('p')
vsT('ap')
axis(1)
dev.off()
#-----------------------------------------------------------------------------#

### Plot a-e
pdf(paste(dir,'a-e.pdf',sep=''),width=10,height=10)
plot(aei[1,,'a'], aei[1,,'e'], 
	xlim=c( 0., 10. ), 
	ylim=c( 0., 1.), 
	pch=20,col=fcol[1])

for (i in 2:length(files))	{
	if (fcol[i] == 'gray')	{
		a = aei[i,-1,'a']
		e = aei[i,-1,'e']
		points( a[e<=1.], e[e<=1.],pch=20,cex=.2, col=fcol[i])
	}}
for (i in 2:length(files))	{
	if (fcol[i] != 'gray')	{
		a = aei[i,-1,'a']
		e = aei[i,-1,'e']
		points( a[e<=1.], e[e<=1.],pch=20,cex=.2, col=fcol[i])
	}}
dev.off()

#-----------------------------------------------------------------------------#
#Function to make barplot of a parameter showing fates of objs
fatebars = function(xcol)	{
	xmin = min(aei[ !isstar, 1, xcol ], na.rm=T)
	xmax = max(aei[ !isstar, 1, xcol ], na.rm=T)
	nbr  = 30
	brk  = xmin + (0:nbr)*((xmax-xmin)/nbr)

	h=array( dim=c( length(allfates), nbr ),dimnames=list(allfates,NULL) )
	for (j in 1:length(allfates))	{
		h[j,]=hist( aei[ ((ult==allfates[j]) & !isstar), 1, xcol ],
			br=brk, plot=F)$counts
		}
	barplot(h, space=0, xlab=xcol, legend.text=T,args.legend=list(x='topleft') )
	axvals = axisTicks( c(xmin,xmax), log=F)
	axlocs = (nbr+1)*(axvals-xmin)/(xmax-xmin)
	axis(side=1,at=axlocs,labels=axvals)
	}
fatebars('a')
#-----------------------------------------------------------------------------#
### Where mass from initial disk ends up

allfates=levels(as.factor(ult))
plotprms = c('a','mass')

n=ceiling(sqrt(length(plotprms)))
if (length(plotprms) <= n*(n-1)) m=n-1 else m=n

pdf(paste(dir,'massfates.pdf',sep=''),height=m*3,width=n*3)
par(mfrow=c(m,n), mar=c(5.1,2.1,1.1,1.1))

for (i in 1:length(plotprms))	{
	fatebars(plotprms[i])
	}
	
dev.off()
##-----------------------------------------------------------------------------#

## stacked barplots of time
#brk=0:10*10^6
#stacktime=function(dest1,br)	{
#	bla=array(dim=c(2,10))
#	brk=0:10*10^6
#	bla[1,]=hist(earthrocks$Time[earthrocks$Destination==dest1],breaks=brk,plot=F)$counts/dim(earthrocks)[1]
#	bla[2,]=hist( marsrocks$Time[ marsrocks$Destination==dest1],breaks=brk,plot=F)$counts/dim( marsrocks)[1]
#	return(bla)	}

#pdf('stackedtime.pdf', height=13, width=4.5)
#par(mfrow=c(6,1), mar=c(2.5,2.5,1.5,0.5), oma=c(2,2,2,0))

#over=c(0,2.75,0,0,9,9,3,3,3)
#up= c(0,.003,0,0,.30,.10,.00012,.00003,.012)

#for (j in c(2,5:9))	{
#thing=stacktime(destlevels[j])
#barplot(thing, col=c('blue','red'), space=0, log='', offset=0,	#, ylim=c(min(thing[thing>0]), max(thing))
#legend.text=c(paste('Earth to',destlevels[j]), paste('Mars to',destlevels[j])), args.legend=list(x=over[j], y=up[j],bg='white'))
#axis(side=1,at=c(0,5,10),labels=c(0,5,10))
#}

#mtext('Time (Myr)', side=1, outer=T)
#mtext('Relative frequency', side=2, outer=T, line=0.5)
#mtext('Transfer over time', side=3, outer=T, line=0.5)
#dev.off()


