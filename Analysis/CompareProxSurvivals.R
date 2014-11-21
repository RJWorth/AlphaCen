

n=34
cases = c('B-2','B-3')
selection = 'B'

### Make array of directory numbers
dirs=array(data=NA,dim = c(n))
for (i in 1:n) dirs[i] = sprintf("%02d", i)

### Get number of cases from dimensions of SurvGrid.txt file
### Note: file also contains 3 footer lines and one column names line

### 3D array for data to be filled into
data = array(data = NA, dim = c(n,length(cases),2), 
	dimnames = list(dirs,cases,c('surv','stab')) )

### Read data
for (i in 1:n)	{
	print(dirs[i])
	GridFile = paste('Proxlike/Prx',dirs[i],'/SurvGrid.txt',sep='')
	GridLen = length(readLines(GridFile))-4
	readin = read.table(GridFile,header=TRUE,row.names=1,nrows=GridLen)
	wantrows = grep(selection,rownames(readin))
	for (j in 1:length(wantrows))	{
		for (k in 1:2)	{
			data[i,j,k] = readin[wantrows[j],k]
		}	# k, surv/stab
	}	# j, disk num.
	}	# i, dirs

pdf('Proxlike/PrxDisksSurvival.pdf')
plot(1:n, data[,1,1], pch=1,col='blue',
	ylim=c(0,1), main='Final Disk Survival Fractions',
	xlab='Prx Simulation #', ylab='Fraction surviving/stable')
	# plot stability in closed circles (20), survival in open (1)
	# plot binaries in blue, triples in red
	for (i in 1:length(cases))	{
		if (length(grep('2',colnames(data)[i]))>0) c='blue' else c='red'
		points(1:n, data[,i,1], pch= 1, col=c)
		points(1:n, data[,i,2], pch=20, col=c)
	}

legend('topright',legend=c('Binary','Triple','Survival','Stability'),
	pch=c(20,20,1,20),
	col=c('blue','red','black','black'))

dev.off()

### Make array of just the average binary and triple stability rates for each
avgs = array(data=NA, dim = c(n, 4), dimnames = list(dirs,c('Bin','Tri','delta','match')))
for (i in 1:length(dirs))	{
	avgs[i,1] = mean(data[i, grep('2',colnames(data)) ,2])
	avgs[i,2] = mean(data[i, grep('3',colnames(data)) ,2])
	avgs[i,3] = mean(avgs[i,1] - avgs[i,2])
	}	# i, dirs

### Read in matching parameter, calculated from Compare.py
### ( := how many of aei parameters match (0-3))
FinalParams = rep(0,length(dirs))
	for (i in 1:length(dirs))	{
		fileName = paste('Proxlike/Prx',dirs[i],'/match.txt',sep='')
		FinalParams[i] = as.integer(readChar(fileName, file.info(fileName)$size))
	}	# i, dirs
avgs[,4]=FinalParams

sink('Proxlike/DiskSummary.txt')
print(avgs)
sink()


