### Calculate % of disk surviving, as a f'n of time
outerbound = min(star[2,,'a'],na.rm=T)
surviving = matrix(, nrow = length(time), ncol = dim(disk)[1])
stable = matrix(, nrow = length(time), ncol = dim(disk)[1])
surv.per = rep(0., length(time))
stab.per = rep(0., length(time))
for (i in 1:length(disknames)) {
	surviving[,i] =  !is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
	stable[,i]    = (!is.na(disk[i,,'r']) & !is.na(disk[i,,'a'])
					 & (abs(disk[i,,'r']-disk[i,1,'r']) <= 1.)
					 & ((cent==2) | (disk[i,,'a'] >= 0.)))
	}
for (j in 1:length(time))	{
	surv.per[j] = sum(surviving[j,])/n
	stab.per[j] = sum(   stable[j,])/n
	}
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
r0=disk[,1,'r']
AxisFraction = (0:length(r0))/length(r0)
LabelSuggestions = 2*(0:10)
AxLabs = LabelSuggestions[LabelSuggestions>min(r0) & LabelSuggestions<max(r0)]
AxLocs = rep(NA,length(AxLabs))
for (i in 1:length(AxLabs))	AxLocs[i] = AxisFraction[FindClosest(r0,AxLabs[i])]

FindClosest(r0,5)



