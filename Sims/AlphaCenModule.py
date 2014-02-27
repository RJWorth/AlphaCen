###############################################################################
### Function to count the number of lines in a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

###############################################################################
### Pick random parameters for the stars and make a new big.in
def MakeBigRand(whichdir,whichtime, cent,
	aBmin, aBmax, eBmin, eBmax, iBmin, iBmax,
	aCmin, aCmax, eCmin, eCmax, iCmin, iCmax):
### Needed modules
	import numpy, os
	from random import random, uniform
	from math import pi, sin, cos
#	from decimal import Decimal, getcontext
#	getcontext().prec = 6

### Constants/variables
	AU   = 1.496e13			# cm/AU
	day  = 24.*3600.		# s/day
	MSun = 1.989e33			# g
	MB   = 0.934			# in MSun
	MC   = 0.123			# in MSun
	rad  = 2*pi/360.		# multiply degrees by this to get radians

	print('MakeBigRand '+whichdir+'/In/big.in, '+str(whichtime))

### Pick random a, e, i for B and C
#	aB = Decimal(uniform(aBmin, aBmax))
	aB = uniform(aBmin, aBmax)
	eB = uniform(eBmin, eBmax)
	iB = uniform(iBmin, iBmax)
	gB = uniform(0.0, 360.0)	#these had been converted to radians as written
	nB = uniform(0.0, 360.0)	# not sure why? program wants degrees
	mB = uniform(0.0, 360.0)

	eC = uniform(eCmin, eCmax)
	iC = uniform(iCmin, iCmax)
	gC = uniform(0.0, 360.0)
	nC = uniform(0.0, 360.0)
	mC = uniform(0.0, 360.0)

	aCfactor=aB*(1+eB)/(1-eC)
	aC = uniform(aCmin*aCfactor, aCmax*aCfactor)

### Write big file
	bigfile=open(whichdir+'/In/big.in','w')
	bigfile.write(")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n")
	bigfile.write(") Lines beginning with `)' are ignored.\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	bigfile.write("style (Cartesian, Asteroidal, Cometary) = Asteroidal\n")
	bigfile.write(" epoch (in days) = 0.0\n")
	bigfile.write(")---------------------------------------------------------------------\n")
	if (cent == 'A'):
		bigfile.write('AlCenB      m=0.934\n')
	if (cent == 'B'):
		bigfile.write('AlCenA      m=1.105\n')
	bigfile.write("  "+repr(aB)+"  "+repr(eB)+"  "+repr(iB)+"\n")
	bigfile.write("  "+repr(gB)+"  "+repr(nB)+"  "+repr(mB)+"\n")
	bigfile.write("  0.0  0.0  0.0\n")

	bigfile.write('PrxCen     m=0.123\n')
	bigfile.write("  "+repr(aC)+"  "+repr(eC)+"  "+repr(iC)+"\n")
	bigfile.write("  "+repr(gC)+"  "+repr(nC)+"  "+repr(mC)+"\n")
	bigfile.write("  0.0  0.0  0.0\n")

	bigfile.close()

	InParams=open(whichdir+'/InParams.txt','a')
	if os.path.getsize(whichdir+'/InParams.txt')==0:
		InParams.write('                 aB                  eB                  iB                  aC                  eC                  iC                 gB                  nB                  mB                  gC                  nC                  mC\n')
	InParams.write(" ".join([repr(aB).rjust(19), repr(eB).rjust(19), repr(iB).rjust(19), repr(aC).rjust(19), repr(eC).rjust(19), repr(iC).rjust(19),
	repr(gB).rjust(19), repr(nB).rjust(19), repr(mB).rjust(19), repr(gC).rjust(19), repr(nC).rjust(19), repr(mC).rjust(19),"\n"]))
	InParams.close()


###############################################################################
# A program to read the important bits and record in summary.txt
def Summary(whichdir,whichtime,cent,tmax,ThisT):
	import numpy, os, AlphaCenModule
	from math import log10
	from operator import add
		
	print('Summary '+whichdir+'/Out, '+whichtime)

### Get input parameters (aei for B and C)
	InParams=open(whichdir+'/InParams.txt','r')
	InParams.seek(-300,2)
	In2=InParams.readlines()[-1]
	In=In2.split()[0:6]
	In1=In2.split()[6:12]
	rnd =[2,4,2]
	rnd6=[1,3,1, 1,3,1]
	space=map(add, [5,3,5, 6,3,5], rnd6)
	In=[str(round(float(In[i]),rnd6[i])).rjust(space[i]) for i in range(6)]
	InParams.close()

### Get output parameters (aei for B and C) from AEI files
 	element=open(whichdir+'/Out/element.out','r')
	nelem=AlphaCenModule.file_len(whichdir+'/Out/element.out')
	elem=element.readlines()
	element.close()

	jstB=map(add, [5,7,6],rnd)
	jstC=map(add, [9,7,6],rnd)
	CnB = ['-'.rjust(jstB[i]) for i in range(3)]
	Prx = ['-'.rjust(jstC[i]) for i in range(3)]
	if (nelem > 5):
		if (elem[5].split()[0] == 'AlCenB'):
			CnB = elem[5].split()[1:4]
			CnB=[str(round(float(CnB[i]),rnd[i])).rjust(jstB[i]) for i in range(3)]
		if (elem[5].split()[0] == 'PrxCen'):
			Prx = elem[5].split()[1:4]
			Prx=[str(round(float(Prx[i]),rnd[i])).rjust(jstC[i]) for i in range(3)]
	if (nelem > 6):
		if (elem[6].split()[0] == 'AlCenB'):
			CnB = elem[6].split()[1:4]
			CnB=[str(round(float(CnB[i]),rnd[i])).rjust(jstB[i]) for i in range(3)]
		if (elem[6].split()[0] == 'PrxCen'):
			Prx = elem[6].split()[1:4]
			Prx=[str(round(float(Prx[i]),rnd[i])).rjust(jstC[i]) for i in range(3)]

### Get last timestep and object's fate from info.out
	name,dest,time = AlphaCenModule.ReadInfo(whichdir)
	DestB,TimeB = '-'.rjust(8),str(ThisT).rjust(13)
	DestC,TimeC = '-'.rjust(8),str(ThisT).rjust(13)

	for j in range(len(name)):
		if (name[j] == 'AlCenB'):
			DestB = dest[j].rjust(8)
			TimeB = time[j].rjust(13)
		if (name[j] == 'PrxCen'):
			DestC = dest[j].rjust(8)
			TimeC = time[j].rjust(13)

### Arrange nicely
	summary=In+CnB+Prx+[str(round(float(TimeB))).rjust(13),DestB,
	str(round(float(TimeC))).rjust(13),DestC,'\n']

### Write to summary.out
### but only if the simulation is ending

### Determine whether an object is missing and the run can be ended
	Btime = (float(TimeB)==tmax)	# has either survived the max sim length?
	Ctime = (float(TimeC)==tmax)	# or
	Bdest = (DestB!='-'.rjust(8))	# has either been ejected etc?
	Cdest = (DestC!='-'.rjust(8))	# (should be mutually exclusive)

	stop=(Btime | Ctime | Bdest | Cdest)
	print('1e'+str(int(log10(ThisT)))+', stop='+
	str(Btime)[0]+str(Ctime)[0]+str(Bdest)[0]+str(Cdest)[0])

	StopFile=open(whichdir+'/stopfile.txt','w')
	StopFile.write(str(stop)+'\n')
	StopFile.close()

### If a Proxima-like C was created, stop the series of runs
	bigstop=False
	if (Btime & Ctime):
		bigstop=float(Prx[0])>5000.	#test=0, usually 5000	
		print('whichtime='+str(whichtime)+', bigstop='+str(bigstop))
	
	BigStopFile=open(whichdir+'/bigstopfile.txt','w')
	BigStopFile.write(str(bigstop)+'\n')
	BigStopFile.close()

### write summary
	sumpath=whichdir+'/summary.out'
	if 'Good' in whichdir:
		sumpath='Good/summary.out'		

	if (stop==True | bigstop==True):
		SumFile=open(sumpath, 'a')
		if os.path.getsize(sumpath)==0:
			SumFile.write('    aB     eB     iB      aC     eC     iC     '+\
		'aB2         eB2      iB2         aC2         eC2      iC2          '+\
		'  tB    destB            tC    destC\n')
		SumFile.write(" ".join(summary))
		SumFile.close()

############################################################################
# A program to read the collision info from info.out
def ReadInfo(whichdir):
	import numpy, os, AlphaCenModule	

 	InfoFile=open(whichdir+'/Out/info.out','r')
	InfoLen=AlphaCenModule.file_len(whichdir+'/Out/info.out')
	AllInfo=InfoFile.readlines()
	InfoFile.close()

	start=where(AllInfo,"   Beginning the main integration.\n")
	ends=where(AllInfo,"   Integration complete.\n")
	nloops=len(ends)

	# Read the last loop
	if (nloops == 1):	
		want=[AllInfo[j] for j in range(start[0]+2,ends[0]-1)]
	else:
		want=[AllInfo[j] for j in range(ends[nloops-2]+10,ends[nloops-1]-1)]
	
	# Extract name, dest, and time from 'want' section
	name,dest,time = [], [], []
	if (len(want)>0):
		name,time,dest = ['' for i in range(len(want))], [
		'' for i in range(len(want))], ['' for i in range(len(want))]
		for j in range(len(want)):
			splitline=want[j].split()
			if len(splitline)==8 and splitline[0]!='Continuing':
				name[j],dest[j],time[j]=splitline[4],splitline[0],splitline[6]
			elif len(splitline)==9:
				name[j],dest[j],time[j]=splitline[0], 'Center',splitline[7]
			elif len(splitline)==5 and splitline[0]!='Fractional':
				name[j],dest[j],time[j]=splitline[0],'ejected',splitline[3]

	return name,dest,time

###############################################################################
# Read data from all directories and put in one file

def SumAll(whichdirs,cent):
	import numpy, os, AlphaCenModule
#	from random import random
#	from math import pi, sin, cos
	
	print('SumAll SumAll.out: '+", ".join(whichdirs))

### Read in summary.out from each directory, compile into SumAll.out
	Sum=[]
	for j in range(len(whichdirs)):
 		DirSumFile=open(whichdirs[j]+'/summary.out','r')
		DirSumLen=AlphaCenModule.file_len(whichdirs[j]+'/summary.out')

		DirSum=DirSumFile.readlines()
		DirSum[0]='Dir '+DirSum[0]
		for k in range(1,len(DirSum)):
			DirSum[k]=str(j+1).rjust(3)+' '+DirSum[k]
		if (len(Sum) < 1):
			Sum=DirSum
		else:
			Sum=Sum+DirSum[1:]
		DirSumFile.close()
	
	SumAll=open('SumAll.out','w')
	for j in range(len(Sum)):
#		print(Sum[j])
		SumAll.write(Sum[j])
	SumAll.close()

### Read in InParams.txt from each directory, compile into AllParams.txt
	Par=[]
	for j in range(len(whichdirs)):
 		DirParFile=open(whichdirs[j]+'/InParams.txt','r')
		DirParLen=AlphaCenModule.file_len(whichdirs[j]+'/InParams.txt')

		DirPar=DirParFile.readlines()
		DirPar[0]='Dir '+DirPar[0]
		for k in range(1,len(DirPar)):
			DirPar[k]=str(j+1).rjust(3)+' '+DirPar[k]
		if (len(Par) < 1):
			Par=DirPar
		else:
			Par=Par+DirPar[1:]
		DirParFile.close()
	
	ParAll=open('AllParams.txt','w')
	for j in range(len(Par)):
		ParAll.write(Par[j])
	ParAll.close()

###############################################################################
# Returns indices where AList==AnElement
def where(AList,AnElement):
	inds=[]
	for j in range(len(AList)):
		if (AList[j]==AnElement):
			inds=inds+[j]
	
	return inds





