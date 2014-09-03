
sums=list()
ncols=array()

###############################################################################
for (i in 1:length(sumfiles))	{
	sums[[i]]=read.table(sumfiles[i], header=T, na.strings='-')
	ncols[i]=dim(sums[[i]])[2]
	}
###############################################################################
SumAll=sums[[1]][,1:min(ncols)]
if (length(sums)>1)	{
	for (i in 2:(length(sums))) SumAll=rbind( SumAll,sums[[i]][,1:min(ncols)] )	

	}	#if
SumAll$tB=10^SumAll$logtB
SumAll$tC=10^SumAll$logtC

###############################################################################

attach(SumAll)

