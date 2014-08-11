### Read the custom aei file (Dir/Out/AeiOutFiles/TimeData.txt) and make plots
# first must create file by running:
# python -c 'import AlphaCenModule as AC; AC.WriteAEI("Dir",1e9)'

### From script
args <- commandArgs(trailingOnly = F)
myarg <- sub("-","",args[length(args)])
Dir=myarg		# which directory to use, as string
### or manual input
Dir='Proxlike/080614/Prx01/Original'

### Read in data
aei=read.table(paste(Dir,'/Out/AeiOutFiles/TimeData.txt',sep=''), header=T)	
	mSun = 2e30
	mA   = 1.105*mSun
	mB   = 0.934*mSun
	mC   = 0.123*mSun
m2 = mA+mB
m3 = m2+mC

EB  = aei$epsB/m2
EC  = aei$epsC/m3
aei = cbind(aei,EB,EC)

attach(aei)

### Make a-e plot
pdf(paste(Dir,'/Out/AeiOutFiles/EvsA.pdf',sep=''), width=7, height=3.5)
par(mfrow=c(1,2), mar=c(4, 4, 1, 1))
plot( aB,eB, type='l', col='orange')
plot( aC,eC, type='l', col='red',    log='')
dev.off()

### Make plots over time of each eps, a, e, and i
pdf(paste(Dir,'/Out/AeiOutFiles/EpsAEIvsT.pdf',sep=''))
par(mfrow=c(2,2), mar=c(4, 4, 1, 1))

#plot(t,epsB-epsB[1], type='l', col='orange', log='x', 
#	ylim=c(min( c(epsB-epsB[1],epsC-epsC[1]) ),
#		   max( c(epsB-epsB[1],epsC-epsC[1]) ) ),
#	xlab='Time (yrs)', ylab='Specific orbital energy change (J/kg?)')
#lines(t,epsC-epsC[1],type='l', col='red')
#abline(h=-epsC[1], col='red',lty=2)

plot(t,EB-EB[1], type='l', col='orange', log='x', 
	ylim=c(min( c(EB-EB[1],EC-EC[1]) ),
		   max( c(EB-EB[1],EC-EC[1]) ) ),
	xlab='Time (yrs)', ylab='Orbital energy change (J?)')
lines(t,EC-EC[1],type='l', col='red')
abline(h=-EC[1], col='red',lty=2)

plot(t,  abs(aC), type='l', col='red', log='xy', ylim=c(min(aB),max(aC)),
	xlab='Time (yrs)', ylab='|Semimajor axis| (AU)')
lines(t, abs(aB), type='l', col='orange')

plot(t,  eC, type='l', col='red', log='x', ylim=c(0,1.5),
	xlab='Time (yrs)', ylab='Eccentricity')
lines(t, eB, type='l', col='orange')
abline(h=1.0, col='red',lty=2)

plot(t,  iC, type='l', col='red', log='x', ylim=c(0,180),
	xlab='Time (yrs)', ylab='Inclination (degrees)')
lines(t, iB, type='l', col='orange')
dev.off()

### Close out, restore defaults
detach(aei)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))

