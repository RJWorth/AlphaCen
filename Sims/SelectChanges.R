
par(mfrow=c(4,3))

meds=array(dim=c(6,2))
	colnames(meds)=c('Early','Late')
	rownames(meds)=c('aB','eB','iB','aC','eC','iC')

for (i in 1:2)	{
ms=c('Early','Late')[i]
dels=list(daB[Case==ms & prox==T], deB[Case==ms & prox==T], diB[Case==ms & prox==T], 
          daC[Case==ms & prox==T], deC[Case==ms & prox==T], diC[Case==ms & prox==T])

for (j in 1:6)	{
	meds[j,i]=median(dels[[j]])
	print(summary(dels[[j]]))
	hist(dels[[j]],br=10)
	}
}

print(meds)

final=c( 23.7 ,0.5179,0, 10000,.9,0)
initial=-meds+cbind(final,final)
print(initial)

