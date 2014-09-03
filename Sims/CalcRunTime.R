###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')

if(dim(runtime)[2]==6) {
	colnames(runtime)=c('Dir','machine','n','user','vers','t')
	nlist=strsplit(as.vector(runtime$n),'/')
	for (i in 1:length(nlist)) if (length(nlist[[i]])<2) nlist[[i]]=c(nlist[[i]],NA)
	runtime$nmax= as.numeric(sapply(nlist, "[[", 2))
	runtime$n   = as.numeric(sapply(nlist, "[[", 1))	}
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
avg  =rep(0.,length(machines))
batch=rep(0.,length(machines))
	for (i in 1:length(machines))	{
	avg[i]=sum(t[machine==machines[i]])/sum(n[machine==machines[i]])
	batch[i]=mean(t[machine==machines[i]])
		}

cat('Runtimes in hrs:\n')
cat(sprintf("%7s",machines))
cat('\n')
cat(sprintf("%7.2f",avg))
cat('\n')
cat('Runtime per batch (days):\n')
cat(sprintf("%7.2f",batch/24.))
cat('\n')

palette(rainbow(length(machines)))

pdf(paste(prefix,'runtimes.pdf',sep=''),height=4,width=8)
par(mfrow=c(1,2))

plot(t/n, col=as.numeric(machine), pch=20, 
	xlab='batch', ylab='Time per sim (hrs)')
for (i in 1:length(machines)) abline(h=avg[i],col=i)
legend('topleft',pch=20,col=1:length(machines), legend=c(machines))

plot(t/n, col=as.numeric(machine), pch=20, 
	xlab='batch', ylab='Time per batch (days)')
for (i in 1:length(machines)) abline(h=batch[i]/24.,col=i)
legend('topleft',pch=20,col=1:length(machines), legend=c(machines))


dev.off()

detach(runtime)

