### Get run version from input variable
args <- commandArgs(trailingOnly = F)
        print(length(args))
        print(args)
#if (length(args) >2) options(echo = FALSE)
basedir <- sub("-","",args[length(args)])
n       =  sub("-","",args[length(args)-1])

### If version not given at input, try using version specified here
if (length(dir) == 0 |
        basedir=="-no-readline" |
        basedir=="/usr/lib/R/bin/exec/R" |
        basedir=="/usr/lib64/R/bin/exec/R" |
        basedir=="/Library/Frameworks/R.framework/Resources/bin/exec/x86_64/R") {
                print('no args')
                basedir='Proxlike/'
		n=34
                }
        print(basedir)
	print(n)
#basedir='Proxlike/'
#n=34
pre='Prx'
cases = c('B-2','B-3')
selection = 'B'

### Make array of directory numbers
dirs=array(data=NA,dim = c(n))
if (basedir=='Proxlike/071714/')	{
	for (i in 1:n) dirs[i] = as.character(i)
} else {
	for (i in 1:n) dirs[i] = sprintf("%02d", i)
}
print(dirs)

### Get number of cases from dimensions of SurvGrid.txt file
### Note: file also contains 3 footer lines and one column names line

### 3D array for data to be filled into
cols=c('surv','stab','r','r/a')
data = array(data = NA, dim = c(n,length(cases),length(cols)), 
	dimnames = list(dirs,cases,cols) )

### Read data
for (i in 1:n)	{
	print(dirs[i])
	GridFile = paste(basedir,pre,dirs[i],'/SurvGrid.txt',sep='')
	GridLen = length(readLines(GridFile))-4
	readin = read.table(GridFile,header=TRUE,row.names=1,nrows=GridLen)
	wantrows = grep(selection,rownames(readin))
	for (j in 1:length(wantrows))	{
		for (k in 1:length(cols))	{
			data[i,j,k] = readin[wantrows[j],k]
		}	# k, surv/stab
	}	# j, disk num.
	}	# i, dirs

### Made index to reorder sims
ind = order(data[,length(cases),length(cols)])	# sort by B-3 stability

pdf(paste(basedir,'PrxDisksSurvival.pdf',sep=''))
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
cols = c('Bin','Tri','delta','match','prox','t','aB','eB','iB','aC','eC','iC')
avgs = array(data=NA, dim = c(n, length(cols)), 
	dimnames = list(dirs[ind],cols))
for (i in 1:length(ind))	{
	avgs[i,1] = mean(data[ind[i], grep('2',colnames(data)) ,2])
	avgs[i,2] = mean(data[ind[i], grep('3',colnames(data)) ,2])
	avgs[i,3] = mean(avgs[i,1] - avgs[i,2])
	}	# i, dirs

### Read in matching parameter and system parameters, calculated from Compare.py
### (matching param := how many of aei parameters match (0-3))
FinalParams = rep(0,length(dirs))
	for (i in 1:length(ind))	{
		fileName = paste(basedir,pre,dirs[ind[i]],'/match.txt',sep='')
		readmatch = readChar(fileName, file.info(fileName)$size)
		words = strsplit(readmatch,' ')[[1]]
		print(readmatch)
		print(words)
		for (j in 1:(length(cols)-3))	{
			print(as.numeric(words[j]))
			avgs[i,j+3] = as.numeric(words[j])
		}
	}	# i, dirs

sink(paste(basedir,'DiskSummary.txt',sep=''))
options(width=100)
print(avgs)
options(width=80)
sink()


