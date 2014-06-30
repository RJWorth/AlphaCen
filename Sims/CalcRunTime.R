###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')
#runtime=read.table(paste(prefix,'runtime.txt',sep=''))
prefix='Plots/'

if(dim(runtime)[2]==6) {
	colnames(runtime)=c('Dir','machine','n','user','vers','t')
	runtime$nmax=as.numeric(sapply(strsplit(as.vector(runtime$n),'/'), "[[", 2))
	runtime$n = as.numeric(sapply(strsplit(as.vector(runtime$n),'/'), "[[", 1))}
if(dim(runtime)[2]==5) colnames(runtime)=c('Dir','machine','n','user','t')
if(dim(runtime)[2]==2) {
	colnames(runtime)=c('Dir','t')
	machine=rep('nova',length(runtime$Dir))
		machine[grep('MDir',runtime$Dir)]='myra'
		machine[grep('CDir',runtime$Dir)]='chloe'
		machine=as.factor(machine)
		n=rep(1,length(runtime$Dir))
		runtime=cbind(runtime,machine,n)	}
runtime$t=runtime$t/3600
attach(runtime)

machines=levels(machine)
avg=rep(0.,length(machines))
	for (i in 1:length(machines))	{
	avg[i]=sum(t[machine==machines[i]])/sum(n[machine==machines[i]])	}

cat('Runtimes in hrs:\n')
cat(sprintf("%7s",machines))
cat('\n')
cat(sprintf("%7.4f",avg))
cat('\n')

palette(rainbow(length(machines)))
pdf(paste(prefix,'runtimes.pdf',sep=''),height=4,width=4)
plot(t/n, col=as.numeric(machine), pch=20, 
	xlab='batch', ylab='Time per sim (hrs)')
for (i in length(machines)) abline(h=avg[i],col=i)
legend('topleft',pch=20,col=1:length(machines), legend=c(machines))
dev.off()

detach(runtime)

