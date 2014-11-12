

n=23	#34
dirs=array(data=NA,dim = c(n))
for (i in 1:n) dirs[i] = sprintf("%02d", i)
data = array(data = NA, dim = c(n,4,2), 
	dimnames = list(dirs,c('A2','A3','B2','B3'),c('surv','stab')) )


for (i in 1:n)	{
	file = paste('Proxlike/Prx',dirs[i],'/SurvGrid.txt',sep='')
	readin = read.table(file,header=TRUE,row.names=1,nrows=4)
	for (j in 1:4)	{
		for (k in 1:2)	{
			data[i,j,k] = readin[j,k]
		}	# k, surv/stab
	}	# j, disk num.
	}	# i, dirs

pdf('Proxlike/PrxDisksSurvival.pdf')
plot(1:n, data[,1,1], pch=1,col='blue',
	ylim=c(0,1), main='Final Disk Survival Fractions',
	xlab='Prx Simulation #', ylab='Fraction surviving/stable')
	points(1:n, data[,3,1], pch=1, col='blue')
	points(1:n, data[,2,1], pch=1, col='red')
	points(1:n, data[,4,1], pch=1, col='red')

	points(1:n, data[,1,2], pch=20, col='blue')
	points(1:n, data[,3,2], pch=20, col='blue')
	points(1:n, data[,2,2], pch=20, col='red')
	points(1:n, data[,4,2], pch=20, col='red')

legend('topright',legend=c('Binary','Triple','Survival','Stability'),
	pch=c(20,20,1,20),
	col=c('blue','red','black','black'))

dev.off()

