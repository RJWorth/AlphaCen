###############################################################################
### Read in wall vs sim time data
times=read.table('looptimes.txt',skip=1)
colnames(times)=c('dir','machine','it','logt','date','time','dt')
attach(times)

###############################################################################
### Make simplified arrays of avg wall time per order of magnitude in sim time
logsimt=3:8
simt=10^logsimt
timevalues <- seq(min(logsimt), max(logsimt), length.out=10)

### Create empty arrays to fill below for each machine
realt=rbind(simt,simt)

struct=list()
exp.model=list()
doubleexp.model=list()
realt.exp1=list()
realt.exp2=list()

### names of machines used
machines=c('chloe','shapiro')

###############################################################################
### For each machine, get average walltime per order of magnitude sim time,
### then fit averages to exponential and double exponential models
### and get times predicted by those models
for (i in 1:length(machines)) {
	for (j in 1:length(simt)) {
	realt[i,j]=mean(dt[ (logt==(j+2)) & (machine==machines[i]) ])
	}
### Reorganize data into structures, for fitting below
struct[[i]] <- structure(list(logsimt, realt[i,]), 
                    .Names = c("logsimt", "realt"), 
                    row.names = c(NA, -length(logsimt)), 
                    class = "data.frame")
### Create models
exp.model[[i]]       <- lm(      log10(realt[i,]) ~ logsimt)
doubleexp.model[[i]] <- lm(log10(log10(realt[i,]))~ logsimt)
### Get predicted data values from models
realt.exp1[[i]] <-     10^(predict(      exp.model[[i]],list(logsimt=timevalues)))
realt.exp2[[i]] <- 10^(10^(predict(doubleexp.model[[i]],list(logsimt=timevalues))))

}	# i

### Make pdf with plot of wall time vs. sim time (shows fit models)
source('Files/PlotFileTimes.R')

detach(times)
