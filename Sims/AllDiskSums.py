
import numpy as np

np.set_printoptions(linewidth=200,precision=3)

suf=['True','07171','072314','081414']

#for i in suf:
i='True'

f = open('Proxlike/Plots/DiskSummary-'+i+'.txt','r')
hdr=np.genfromtxt('Proxlike/Plots/DiskSummary-'+i+'.txt',invalid_raise=False,dtype='S')
hdr=hdr.tolist()
hdr.insert(0,'Num')
thistable = np.genfromtxt('Proxlike/Plots/DiskSummary-'+i+'.txt',
	dtype=['i']+['f' for i in range(len(hdr)-1)],skip_header=1,names=hdr)
	
fname='Proxlike/Plots/DiskSummary-'+i+'.txt'
d=['i']+['f' for j in range(len(hdr)-1)]
np.loadtxt(fname,dtype=d,skiprows=1)
	

