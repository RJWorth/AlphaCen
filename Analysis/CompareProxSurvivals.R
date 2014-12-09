

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

### Made index to reorder sims
ind = order(data[,length(cases),2])	# sort by B-3 stability

pdf('Proxlike/PrxDisksSurvival.pdf')
plot(1:n, data[ind,1,1], pch=1,col='blue', type='n',
	xaxt='n', ylim=c(0,1), 
	main='Final Disk Survival Fractions',
	xlab='Prx Simulation', ylab='Fraction surviving/stable')
	# plot stability in closed circles (20), survival in open (1)
	# plot binaries in blue, triples in red
	for (i in 1:length(cases))	{
		if (length(grep('2',colnames(data)[i]))>0) c='blue' else c='red'
		points(1:n, data[ind,i,1], pch= 1, col=c)
		lines( 1:n, data[ind,i,1], lty= 3, col=c)
		points(1:n, data[ind,i,2], pch=20, col=c)
		lines( 1:n, data[ind,i,2], lty= 1, col=c)
	}

legend('topright',legend=c('Binary','Triple','Survival','Stability'),
	pch=c(20,20,1,20),
	col=c('blue','red','black','black'))

dev.off()

### Make array of just the average binary and triple stability rates for each
avgs = array(data=NA, dim = c(n, 11), 
	dimnames = list(dirs[ind],c('Bin','Tri','delta','match','t','aB','eB','iB','aC','eC','iC')))
for (i in 1:length(ind))	{
	avgs[i,1] = mean(data[ind[i], grep('2',colnames(data)) ,2])
	avgs[i,2] = mean(data[ind[i], grep('3',colnames(data)) ,2])
	avgs[i,3] = mean(avgs[i,1] - avgs[i,2])
	}	# i, dirs

### Read in matching parameter and system parameters, calculated from Compare.py
### (matching param := how many of aei parameters match (0-3))
FinalParams = rep(0,length(dirs))
	for (i in 1:length(ind))	{
		fileName = paste('Proxlike/Prx',dirs[ind[i]],'/match.txt',sep='')
		readmatch = readChar(fileName, file.info(fileName)$size)
		words = strsplit(readmatch,' ')[[1]]
		print(readmatch)
		print(words)
		for (j in 1:8)	{
			print(as.numeric(words[j]))
			avgs[i,j+3] = as.numeric(words[j])
		}
	}	# i, dirs

sink('Proxlike/DiskSummary.txt')
options(width=100)
print(avgs)
options(width=80)
sink()


