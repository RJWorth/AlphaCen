disktimes=read.table('disktimes.txt')
	colnames(disktimes)=c('dir','logt','wallt')
attach(disktimes)

logst=3:7
st=10^logst
wt=rep(NA,length(st))

for (i in 1:length(st)) wt[i] = max(wallt[ logt==logst[i] ])
logwt=log10(wt)

fit = lm(logwt ~ logst)
b = fit$coef[1]
m = fit$coef[2]

### predictions
logpst = 3:9
pst = 10^logpst
pwt = 10^b * pst^m

print(paste('1e6 sims done in',pwt[4]/3600/24,'days'))
print(paste('1e7 sims done in',pwt[5]/3600/24,'days'))
print(paste('1e8 sims done in',pwt[6]/3600/24,'days'))
print(paste('1e9 sims done in',pwt[7]/3600/24,'days'))

disktimes2 = read.table('disktimes2.txt')
	colnames(disktimes2)=c('dir2','logt2','wallt2')
attach(disktimes2)

plot(10^logt,wallt,pch=20,log='xy')
abline(fit,col='red')
points(10^logt2,wallt2,pch=20,col='blue')


