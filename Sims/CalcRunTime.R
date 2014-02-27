###############################################################################
### Calculate average runtime per sim
runtime=read.table('runtime.txt')
#runtime=read.table(paste(prefix,'runtime.txt',sep=''))

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

avg=c(0,0,0)
	avg[1]=sum(t[machine=='nova'])/sum(n[machine=='nova'])
	avg[2]=sum(t[machine=='myra'])/sum(n[machine=='myra'])
	avg[3]=sum(t[machine=='chloe'])/sum(n[machine=='chloe'])

print('runtimes in hrs')
print('   n         m         c')
print(avg)

palette(c('blue','green','red'))
pdf(paste(prefix,'runtimes.pdf',sep=''),height=4,width=4)
plot(t/n,col=machine,pch=20, xlab='batch',ylab='Time per sim (hrs)')
abline(h=avg[1],col=3)
abline(h=avg[2],col=2)
abline(h=avg[3],col=1)
legend('topleft',pch=20,col=c(3,2,1), legend=c('nova','myra','chloe'))
dev.off()

detach(runtime)

