
print('Analyzing data')
### Radius key
r0=disk[,1,'r']
aBin = star[2,1,'a']	# initial binary separation

### Calculate % of disk surviving and location of outer boundary,
### as f'ns of time
outerbound = min(star[2,,'a'],na.rm=T)
surviving = matrix(, nrow = length(time), ncol = dim(disk)[1])
stable    = matrix(, nrow = length(time), ncol = dim(disk)[1])
above     = matrix(, nrow = length(time), ncol = dim(disk)[1])
below     = matrix(, nrow = length(time), ncol = dim(disk)[1])
side.fits = matrix(, nrow = length(time), ncol = dim(disk)[1])
side.smth = matrix(, nrow = length(time), ncol = dim(disk)[1])
surv.per = rep(0., length(time))
stab.per = rep(0., length(time))
edge     = rep(0., length(time))
# calculate whether each obj is surviving/stable at each timestep
for (i in 1:length(disknames)) {
	surviving[,i] =  !is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
	stable[,i]    = (!is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
					 & (abs(disk[i,,'r']-disk[i,1,'r']) <= 1.)
					 & ((cent==2) | (disk[i,,'a'] >= 0.)))
	} #i
# Calculate avg stability on either side of each obj at each timestep
for (i in 1:length(disknames)) {
for (j in 1:length(time))	{
	below[j,i] = mean(stable[j,1:i])
	if (i < n)	above[j,i] = mean(stable[j,(i+1):n]) else {
		if (stable[j,i]==TRUE) above[j,i]=1 else above[j,i]=0 }
	} #j
	} #i
# how closely does each point fit a solid disk with nothing beyond?
side.fits = (1-below) + (above)
# smooth fit by taking 3-point moving average
for (j in 1:length(time)) side.smth[j,] = MovingAvg(side.fits[j,],1)
# Calculate percent of objs surviving at each step; fit edge of surviving disk
for (j in 1:length(time))	{
	surv.per[j] = sum(surviving[j,])/n
	stab.per[j] = sum(   stable[j,])/n
	if ( all( side.fits[j,]==1 )) edge[j]=100 else{
	edge[j]     = mean( which( side.smth[j,]==min(side.smth[j,],na.rm=T) ) ) }
	} #j
#extent=max(abs(c(disk[,,4:5],star[,,4:5])),na.rm=T)
extent=max(abs(c(disk[,1,4:5],star[,,4:5])),na.rm=T)
grays = gray.colors(n, start = 0., end = 0.8)
heat  = heat.colors(n)

### Resample log scale for time
samprate=100
logt = min(log10(time[-1])) + 
		(max(log10(time))-min(log10(time[-1])))*(0:(samprate-1))/(samprate-1)
resampt = 10^logt
### OR linear scale
#resampt = min(time) + 
#			(max(time)-min(time))*(0:(samprate-1))/(samprate-1)


### Make image matrix of disk survival
diskimg = matrix(, nrow = samprate, ncol = dim(disk)[1])
	for (j in 1:n)	{
		for (k in 1:samprate)	{
			i = which( abs(time-resampt[k]) == min(abs(time-resampt[k])) )
			if (length(i) > 1) i=min(i)
			if ((surviving[i,j]==TRUE) & (stable[i,j] == TRUE)) {
				diskimg[k,j]=j
			} else if ((surviving[i,j]==TRUE) & (stable[i,j] == FALSE)) {
				diskimg[k,j]=1+n
			} else if ((surviving[i,j]==FALSE)) {
				diskimg[k,j]=2+n
			} else {print(paste(k,j))}
	}}

### Measure the difference between actual and expected velocity, 
### as a fraction of expected
#verr=(v(1:100,1)-vorb_expect(1:100,1))/vorb_expect(1:100,1)

### Make a 2D array of initial parameters for each object
#pairframe=cbind(disk[,1,],verr)
pairframe=disk[,1,]

time.all  =               1:maxTlength
timeslice = (maxTlength-99):maxTlength

### Find where to put axis ticks in DiskSurv plot
rA = r0/aBin
AxisFraction = (0:length(rA))/length(rA)
LabelSuggestions = .1*(0:20)
AxLabs = LabelSuggestions[LabelSuggestions>min(rA) & LabelSuggestions<max(rA)]
AxLocs = rep(NA,length(AxLabs))

AxLocs = (approx(rA,(0:(n-1))/(n-1),xout=AxLabs)$y)

#for (i in 1:length(AxLabs))	AxLocs[i] = AxisFraction[FindClosest(rA,AxLabs[i])]


