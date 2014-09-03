###############################################################################
### Make simplified arrays of avg wall time per order of magnitude in sim time
logsimt=3:9
simt=10^logsimt
timevalues <- seq(min(logsimt), max(logsimt), length.out=(length(simt)*2-1))
### Break point between the two models, log(t)
brk=6
### Indices for the short sims before the break, and long sims after
i.shrt = logsimt <= brk
i.long = logsimt >= brk

t.shrt = timevalues <= brk
t.long = timevalues >= brk

### Read in wall vs sim time data
times=read.table('looptimes1.txt',skip=1)
colnames(times)=c('dir','machine','it','logt','date','time','dt')
attach(times)

new=read.table('looptimes.txt',skip=1)
colnames(new)=c('dir','machine','it','logt','date','time','b','dt')

### Indices for rows of times that are before/after break
r.shrt = logt <= brk
r.long = logt >= brk
times1=times[r.shrt,]
times2=times[r.long,]

# For shorter times (logt<=5) a double exponential fits better, but 
# for logt>5, a regular exponential is best

###############################################################################

### Create empty arrays to fill below for each machine
realt=rbind(simt,simt)

shrtsims=list()
longsims=list()

lin.model=list()
exp.model=list()
dbl.model=list()

realt.lin=list()
realt.exp=list()
realt.dbl=list()

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
shrtsims[[i]] <- structure(list(logsimt[i.shrt], realt[i,i.shrt]), 
                    .Names = c("logsimt", "realt"), 
                    row.names = c(NA, -sum(i.shrt)), 
                    class = "data.frame")
longsims[[i]] <- structure(list(logsimt[i.long], simt[i.long], realt[i,i.long]), 
                    .Names = c("logsimt", "simt", "realt"), 
                    row.names = c(NA, -sum(i.long)), 
                    class = "data.frame")
### Create models
# double-exponential for short sims
dbl.model[[i]] <- lm(log10(log10(realt))~ logsimt, data=shrtsims[[i]],
						weights=c( rep(1,sum(i.shrt)-1), 5) )

# linear and exponential models for long sims
exp.model[[i]] <- lm(      log10(realt) ~ logsimt, data=longsims[[i]] )
lin.model[[i]] <- lm(            realt  ~    simt, data=longsims[[i]] )

### Get predicted data values from models
realt.dbl[[i]] <- 10^(10^(predict(dbl.model[[i]],
									list(logsimt=   timevalues[t.shrt]) ) ) )

realt.exp[[i]] <-     10^(predict(exp.model[[i]],
									list(logsimt=   timevalues[t.long]) ) )
realt.lin[[i]] <-        (predict(lin.model[[i]],
									list(   simt=10^timevalues[t.long]) ) )

}	# i

### Get coefficients from the models in a list
### coef[[machine]][[model]][coefficient]
p=3	#precision for saved coefficient values
coef=list()
	coef[[1]]=list()
	coef[[2]]=list()
c.f=coef
for (i in 1:length(machines)) {
coef[[i]][[1]]=summary(dbl.model[[i]])$coefficients[,1]
	coef[[i]][[1]] = c(10^(10^(coef[[i]][[1]][1])), coef[[i]][[1]][2])
coef[[i]][[2]]=summary(exp.model[[i]])$coefficients[,1]
	coef[[i]][[2]] = c(   (10^(coef[[i]][[2]][1])), coef[[i]][[2]][2])
coef[[i]][[3]]=summary(lin.model[[i]])$coefficients[,1]
	coef[[i]][[3]] = c(   (   (coef[[i]][[3]][1])), coef[[i]][[3]][2])
	}
# formatted coefficients for legend:
for (i in 1:length(machines)){
	for (j in 1:3)	c.f[[i]][[j]]=signif( coef[[i]][[j]], p)	}

### Prediction function based on models
PredictWTfromB = function(bytes, coef.all)	{
	# For chloe, simt>1e6
	C   = coef.all[[1]][[2]][1]
	ex  = coef.all[[1]][[2]][2]
	int = -126072
	m   =   66773

	wt  = C * (10^(-int*ex/m)) * (10^(ex*bytes/m))
	return(wt)	}

PredictWallT.dbl = function(simt, coef)	{
	WallT = coef[1]^(simt^(coef[2]))
	return(WallT)	}
PredictWallT.exp = function(simt, coef)	{
	WallT = coef[1] * simt^(coef[2])
	return(WallT)	}
PredictWallT = function(simt, machine, coef.all)	{
	if (machine == 'shapiro') m=2 else m=1
	WallT = rep(0,length(simt))
	for (i in 1:length(simt))	{
		if (log10(simt[i]) > 6)	{
			WallT[i] = PredictWallT.exp(simt[i], coef.all[[m]][[2]])	
		} else {
			WallT[i] = PredictWallT.dbl(simt[i], coef.all[[m]][[1]])
#			WallT.dbl = PredictWallT.dbl(simt[i], coef.all[[m]][[1]])
#			WallT.exp = PredictWallT.exp(simt[i], coef.all[[m]][[2]])	
#			WallT[i] = max(WallT.exp, WallT.dbl)
		} # else
	}	# for i
	return(WallT)
	}	# fn PredictWallT
n=100.
faket=10^(((3*n):(9*n))/n)

fakewt.c=PredictWallT(faket,  'chloe',coef)
fakewt.s=PredictWallT(faket,'shapiro',coef)

### Write the best fit coefficients to a file
coef.array=rbind(coef[[1]][[1]], coef[[2]][[1]],
				 coef[[1]][[2]], coef[[2]][[2]])
	rownames(coef.array)=c( 'c.dbl','s.dbl','c.exp','s.exp' )
	colnames(coef.array)=c( 'c','exponent' )
sink('RunTimeCoefs.txt')
print(coef.array)
sink()

### Make pdf with plot of wall time vs. sim time (shows fit models)
source('Files/PlotFileTimes.R')

detach(times)
