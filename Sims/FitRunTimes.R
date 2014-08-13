times=read.table('looptimes.txt',skip=1)
colnames(times)=c('dir','machine','it','logt','date','time','dt')
attach(times)

logsimt=3:8
simt=10^logsimt
realt=simt
for (i in 1:length(simt)) realt[i]=mean(dt[logt==(i+2)])

struct <- structure(list(logsimt, realt), 
					.Names = c("logsimt", "realt"), 
					row.names = c(NA, -length(logsimt)), 
					class = "data.frame")


exp.model       <- lm(      log10(realt) ~ logsimt)
doubleexp.model <- lm(log10(log10(realt))~ logsimt)

timevalues <- seq(min(logsimt), max(logsimt), length.out=10)
realt.exp1 <-     10^(predict(      exp.model,list(logsimt=timevalues)))
realt.exp2 <- 10^(10^(predict(doubleexp.model,list(logsimt=timevalues))))

plot(10^logt,dt, pch=20, log='xy', xlab='sim time (yrs)',ylab='walltime (s)')
points(10^logsimt, realt, pch=20,col='red')
lines(10^timevalues, realt.exp1, col = "blue")
lines(10^timevalues, realt.exp2, col = "green")

legend('topleft',legend=c('actual times','average times','walltime ~ (sim t)^n','walltime ~ m^((sim t)^n)'),
	col=c('black','red','blue','green'),lty=c(0,0,1,1),pch=c(20,20,NA,NA))

detach(times)
