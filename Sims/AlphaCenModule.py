###############################################################################
def FileLength(fname):
	'''Function to count the number of lines in a file'''
	import os

	if os.path.getsize(fname)==0.:
		i = -2
	else:
	    with open(fname) as f:
	        for i, l in enumerate(f):
	            pass

	return i + 1

###############################################################################
def where(AList,AnElement):
	'''Returns indices where AList==AnElement'''
	inds=[]
	for j in range(len(AList)):
		if (AList[j]==AnElement):
			inds=inds+[j]
	
	return inds

############################################################################
def WriteObjInFile(WhichDir,names,filename,Header,FirstLines,xv,s):
	'''Write big.in or small.in file'''

	infile=open(WhichDir+'/In/'+filename+'.in','w')
### Header
	for i in range(len(Header)):
		infile.write(Header[i])
### Data
	for i in range(len(names)):
		infile.write(FirstLines[i])
		infile.write("  "+xv[i][0]+"  "+xv[i][1]+"  "+xv[i][2]+"\n")
		infile.write("  "+xv[i][3]+"  "+xv[i][4]+"  "+xv[i][5]+"\n")
		infile.write(s[i])
	infile.close()

###############################################################################
def MakeBigRand(WhichDir,WhichTime, cent,
	aBmin, aBmax, eBmin, eBmax, iBmin, iBmax,
	aCmin, aCmax, eCmin, eCmax, iCmin, iCmax, 
	mB=0.934, mC=0.123):
	'''Pick random parameters for the stars and make a new big.in'''

	print('	MakeBigRand  '+WhichDir+'/In/big.in,                          '+
		  str(WhichTime))

### Needed modules
	import AlphaCenModule as AC
	import numpy, os
	from random import random, uniform
	from math import pi, sin, cos
#	from decimal import Decimal, getcontext
#	getcontext().prec = 6

### Constants
	AU   = 1.496e13			# cm/AU
	day  = 24.*3600.		# s/day
	MSun = 1.989e33			# g
	MB   = mB				# in MSun
	MC   = mC				# in MSun
	rad  = 2*pi/360.		# multiply degrees by this to get radians

### Pick random a, e, i and g, n, m for B and C
	aB = uniform(aBmin, aBmax)
	eB = uniform(eBmin, eBmax)
	iB = uniform(iBmin, iBmax)
	gB, nB, mB = uniform(0.0, 360.0), uniform(0.0, 360.0), uniform(0.0, 360.0)

	eC = uniform(eCmin, eCmax)
	iC = uniform(iCmin, iCmax)
	gC, nC, mC = uniform(0.0, 360.0), uniform(0.0, 360.0), uniform(0.0, 360.0)

	aCfactor = aB*(1+eB)/(1-eC)
	aC = uniform(aCmin*aCfactor, aCmax*aCfactor)

	aei = [[repr(aB), repr(eB), repr(iB), repr(gB), repr(nB), repr(mB)], 
		   [repr(aC), repr(eC), repr(iC), repr(gC), repr(nC), repr(mC)]]

### Read generic big.in file header	
	BigHeadFile=open('BigHeader.txt','r')
	BigHeader=BigHeadFile.readlines()
	BigHeadFile.close()

### First lines for each object
	if (cent == 'A'):
		BigFirstLines=(['AlCenB      m=0.934  r=3.0\n'])
		names=['AlCenB']
	if (cent == 'B'):
		BigFirstLines=(['AlCenA      m=1.105  r=3.0\n'])
		names=['AlCenA']
	BigFirstLines.append('PrxCen     m=0.123  r=3.0\n')
	names.append('PrxCen')

### Spin
### No spin for all objects
	BigS=["  0.0  0.0  0.0\n" for i in range(len(BigFirstLines))]

### Write big file
	AC.WriteObjInFile(
	WhichDir,names,'big',BigHeader,BigFirstLines,aei,BigS)

### Save initial parameters with full precision in InParams.txt
	InParams=open(WhichDir+'/InParams.txt','a')
	if os.path.getsize(WhichDir+'/InParams.txt')==0:
		InParams.write('                 aB                  eB'+\
		'                  iB                  aC                  eC'+\
		'                  iC                 gB                  nB'+\
		'                  mB                  gC                  nC'+\
		'                  mC\n')
	InParams.write(" ".join([repr(aB).rjust(19), repr(eB).rjust(19), 
													repr(iB).rjust(19), 
	repr(aC).rjust(19), repr(eC).rjust(19), repr(iC).rjust(19),
	repr(gB).rjust(19), repr(nB).rjust(19), repr(mB).rjust(19), 
	repr(gC).rjust(19), repr(nC).rjust(19), repr(mC).rjust(19),"\n"]))
	InParams.close()

###############################################################################
def InitParams(WhichDir):
	'''Read initial orbital parameters'''

### Modules needed
	import numpy as np
	from operator import add
	import os 

### Digits
	iRnd =[1,3,1, 1,3,1]	# dec places to round starting a/e/i to, for B/C
	iJst = map(add, [4,2,3, 5,2,4], iRnd) # iRnd+iJst = # char per entry

### Get input parameters (aei for B and C)
	if os.path.getsize(WhichDir+'/InParams.txt')==0:
		print('No InParams file!')
		aeiIn=['-','-','-','-','-','-']
	else:
		InParamsFile=open(WhichDir+'/InParams.txt','r')
		InParamsFile.seek(-300,2)
		RawIn=InParamsFile.readlines()[-1]
		InParamsFile.close()
		aeiIn=RawIn.split()[0:6]	# aei for B, C
#	gnmIn=RawIn.split()[6:12]	# angles for B, C
	# spaces needed for initial parameters = max digits+1+dec points
		aeiIn=[str(round(float(aeiIn[i]),iRnd[i])) for i in range(6)]
	aeiIn=[aeiIn[i].rjust(iJst[i]) for i in range(6)]

	return(aeiIn)

###############################################################################
def Elem(WhichDir):
	'''Get final orbits from element.out'''

### Modules needed
	import numpy as np
	from operator import add
	import AlphaCenModule as AC

### Constants
	fRnd =[2,4,2]			# dec places to round final a/e/i to

### Get output parameters (aei for B and C) from AEI files
 	ElemFile=open(WhichDir+'/Out/element.out','r')
	nElem=AC.FileLength(WhichDir+'/Out/element.out')
	elem=ElemFile.readlines()
	ElemFile.close()
### Check for '*****' entries, trigger bigstop if found
#	bigstop=False
#	if ("*" in elem[5].split()[1:4]):
#		bigstop=True
#		print('element.out contains *** entries')
#		BigStopFile=open(WhichDir+'/bigstopfile.txt','w')
#		BigStopFile.write(str(bigstop)+'\n')
#		BigStopFile.close()

### create vectors to hold a,e,i, x,y,z for CnB and Prx
	jstB=map(add, [5,7,6],fRnd)		# digits+dec points
	jstC=map(add, [9,7,6],fRnd)
	B = ['-'.rjust(jstB[i]) for i in range(3)]
	C = ['-'.rjust(jstC[i]) for i in range(3)]

### if element file contains objects
	if (nElem > 5):
### if the first object is CnB, read its parameters
		if (elem[5].split()[0] == 'AlCenB'):
			B = elem[5].split()[1:4]
			B = [str( round(float(B[i]),fRnd[i]) ).rjust(jstB[i]) 
               for i in range(3)]
### or if the first object is Prx, read its parameters
		if (elem[5].split()[0] == 'PrxCen'):
			C = elem[5].split()[1:4]
			C = [str( round(float(C[i]),fRnd[i]) ).rjust(jstC[i]) 
                 for i in range(3)]
### do the same for the second element
	if (nElem > 6):
		if (elem[6].split()[0] == 'AlCenB'):
			B = elem[6].split()[1:4]
			B = [str( round(float(B[i]),fRnd[i]) ).rjust(jstB[i]) 
                for i in range(3)]
		if (elem[6].split()[0] == 'PrxCen'):
			C = elem[6].split()[1:4]
			C = [str( round(float(C[i]),fRnd[i]) ).rjust(jstC[i]) 
                 for i in range(3)]

	return(B, C)

###############################################################################
def ReadAei(whichdir, filename, index1=-1, index2=0):
	'''Read .aei file to get xyz and uvw (in AU and m/s respectively)'''

	print('	ReadAei      '+whichdir+'/Out/AeiOutFiles/'+filename+', '+\
		   str(index1)+':'+str(index2))

### Modules needed
	import numpy as np

### Get last positions for surviving objects
	aeiFile=open(whichdir+'/Out/AeiOutFiles/'+filename+'.aei','r')
	aei = aeiFile.readlines()
	if index2!=0:
		xv = aei[index1:index2]
	else:
		xv = aei[index1:]
### Arrange xyz uvw data into array
	xv = np.array([i.split()[6:12] for i in xv])

### Convert from strings to floats
	xv = xv.astype(np.float)

	return(xv)

###############################################################################
def GetT(whichdir, filename, index1=-1, index2=0):
	'''Read .aei file to get xyz and uvw (in AU and m/s respectively)'''

	print('	GetT         '+whichdir+'/Out/AeiOutFiles/'+filename+', '+\
		   str(index1)+':'+str(index2))

### Modules needed
	import numpy as np

### Get last positions for surviving objects
	aeiFile=open(whichdir+'/Out/AeiOutFiles/'+filename+'.aei','r')
	aei = aeiFile.readlines()
	if index2!=0:
		xv = aei[index1:index2]
	else:
		xv = aei[index1:]
### Arrange xyz uvw data into array
	t = np.array([i.split()[0] for i in xv])

### Convert from strings to floats
	t = t.astype(np.float)

	return(t)

###############################################################################
def Distance(xv1, xv2):
	'''Get distance between two points'''

	r = (  (xv1[:,0]-xv2[:,0])**2
		 + (xv1[:,1]-xv2[:,1])**2
		 + (xv1[:,2]-xv2[:,2])**2	)**0.5

	return(r)

###############################################################################
def XVtoR(xv):
	'''Take xv vector(s), return r'''

	r = ( xv[:,0]**2 + xv[:,1]**2 + xv[:,2]**2 )**0.5

	return(r)

###############################################################################
def XVtoV(xv):
	'''Take xv vector(s), return v'''

	v = ( xv[:,3]**2 + xv[:,4]**2 + xv[:,5]**2 )**0.5

	return(v)

###############################################################################
def RVtoA(r, v, mu):
	'''Estimate a from v(r) function
	where r is the distance between two stars, v is the relative velocity'''

	a = 1/( 2/r -  v**2/mu )

	return(a)

###############################################################################
def RAtoV(r, a, mu):
	'''Estimate v from v(r) function
	where r is the distance between two stars, v is the relative velocity'''

	v = (mu*( 2/r -  1/a) )**0.5

	return(v)

###############################################################################
def AUtoMKS(xv_AU):
	'''Convert xv in AU, AU/day to mks'''

	import numpy as np
	from mks_constants import AU, day

	xv_mks = np.zeros_like(xv_AU)
	
	if (len(np.shape(xv_AU))==3):
		xv_mks[:,:,0:3] = xv_AU[:,:,0:3]*AU
		xv_mks[:,:,3:6] = xv_AU[:,:,3:6]*AU/day
	elif (len(np.shape(xv_AU))==2):
		xv_mks[:,0:3] = xv_AU[:,0:3]*AU
		xv_mks[:,3:6] = xv_AU[:,3:6]*AU/day

	return(xv_mks)

###############################################################################
def FindCM(m, xv):
	'''Take position and velocity data, 
	give center of mass/momentum versions'''

	from mks_constants import mSun	

	print('	FindCM       m = '+str([('% 1.3g' % (i/mSun)) for i in m]))

### Modules needed
	from operator import add
	import numpy as np

### Constants
	M = sum(m)

### Number of timesteps being looked at
	nobjs =xv.shape[0]
	nsteps=xv.shape[1]

### Create empty array with same shape
	xvCM  = np.zeros_like(xv[0,:,:])

### Calculate transforms
	for j in range(nsteps):
		for k in range(0,6):
			xvCM[j,k] = (1/M)*(sum( m[:]*xv[:,j,k] ))

	return(xvCM)

###############################################################################
def wrtCM(xv, xvCM):
	'''Transform given coordinates to CM frame'''

#	print('	wrtCM')

### Modules needed
	from operator import add
	import numpy as np	

### Constants
#	mSun	= 1.9891e30		# kg
#	m = m/mSun
#	M = sum(m)

### Number of timesteps being looked at
	nobjs =xv.shape[0]
	nsteps=xv.shape[1]

### Create more empty arrays with same shape, to fill later
	xv2CM = np.zeros_like(xv)

### Transform B and C coordinates by adding the transform vectors
	for i in range(nobjs):
		for j in range(nsteps):
			for k in range(6):
				xv2CM[i,j,k] = xv[i,j,k]-xvCM[j,k]

	return(xv2CM)

###############################################################################
def Kinetic(m, v):
	'''Calculate kinetic energy
	where r is the distance between two stars, v is the relative velocity'''

	import numpy as np

	K = [(1./2.)*m[i]*v[i,:]**2. for i in range(len(v[:,0]))]

	return(np.array(K))

###############################################################################
def Potential(m1, m2, r):
	'''Calculate potential energy of one object due to another
	where r is the distance between two stars, v is the relative velocity'''

	from mks_constants import G

	U = -G*m1*m2/r
	U = U/2.		# correction to avoid counting U twice

	return(U)

###############################################################################
def Eps(r, v, mu):
	'''Calculate specific orbital energy
	where r is the distance between two stars, v is the relative velocity.
	Final units should be J/kg or m^2/s^2.'''

	eps = v**2./2.-mu/r

	return(eps)

###############################################################################
def mr(m):
	'''Calculate reduced mass'''

### Needed modules
	import numpy as np

	m_1  = [1./i for i in m]
	mr_1 = sum(m_1)
	mr   = 1./mr_1

	return(mr)

###############################################################################
def h(r, v):
	'''Calculate angular momentum'''

### Needed modules
	import numpy as np

	h = np.cross(r,v)

	return(np.array(h))

###############################################################################
def a(eps, mu):
	'''Calculate semimajor axis of an object's orbit'''

### Needed modules
	import numpy as np

	a = -mu/(2.*eps)

	return(np.array(a))

###############################################################################
#def a(E, m, mu):
#	'''Calculate semimajor axis of an object's orbit'''

#	from mks_constants import G
### Needed modules
#	import numpy as np
#	import AlphaCenModule as AC

#	mr = AC.mr(m)	# reduced mass

#	a = -mu/(2.*(E/mr))

#	return(np.array(a))

###############################################################################
def e(a, eps, hbar, mu):
	'''Calculate eccentricity of an object's orbit'''

### Needed modules
	import numpy as np

#	print(2.*eps*hbar**2)
#	print(mu**2)

#	print(2.*eps*hbar**2/mu**2)
	e = (1.+2.*eps*hbar**2/mu**2)**0.5

	return(np.array(e))

###############################################################################
def i(hz, hbar):
	'''Calculate inclination of an object's orbit (in degrees)'''

### Needed modules
	import numpy as np
	from numpy import arccos, pi

	i = arccos(hz/hbar)*180./pi

	return(np.array(i))

###############################################################################
def GetFinalData(WhichDir,ThisT,mode):
	'''Read in data, calculate derived values, and output for 
	writing/plotting/etc.'''

	print('	GetFinalData '+WhichDir)

### Needed modules
	import numpy as np
#	import os
	import AlphaCenModule as AC
#	from numpy import log10, sqrt, sin, pi
#	from operator import add
#	import subprocess
#	import rvtest as rv
	from mks_constants import G, mSun, AU, day, m, mu
	from numpy import log10

	assert (mode=='triple') | (mode=='binary')

	if (mode == 'binary'):
		m=[m[0], m[1], 0.]
		filenames=['AlCenB']
	elif (mode == 'triple'):
		filenames=['AlCenB', 'PrxCen']

### Number of objects
	nobjs=len(filenames)+1
### Number of timesteps
	ntB=AC.FileLength(WhichDir+'/Out/AeiOutFiles/'+filenames[0]+'.aei')-4
	if (mode == 'triple'):
		ntC=AC.FileLength(WhichDir+'/Out/AeiOutFiles/'+filenames[1]+'.aei')-4
	else:
		ntC=0.
	
	t = AC.GetT(WhichDir, filenames[0], 4, ntB+4)

	Bind = [4,ntB+4]
	Cind = [4,ntC+4]

########## Get object's fate and collision/ejection time from info.out #######
	name,dest,time = AC.ReadInfo(WhichDir)
	DestB,TimeB = '-'.rjust(8),str(ThisT).rjust(13)
	DestC,TimeC = '-'.rjust(8),str(ThisT).rjust(13)

	for j in range(len(name)):
		if (name[j] == 'AlCenB'):
			DestB = dest[j].rjust(8)
			TimeB = time[j].rjust(13)
		if (name[j] == 'PrxCen'):
			DestC = dest[j].rjust(8)
			TimeC = time[j].rjust(13)

############# Read in original x and v values ###########################
### NOTATION:
# 3D array: xvA = xv of all objects w.r.t. CMA (in AU units)
# xvA[i,j,k] = xv of object i, time j, column k
# Units are mks unless specified otherwise by tacking _AU on the end, 
# then it's in the AU units that MERCURY outputs (x=AU, v=AU/day).

# rAB = distance between A and B
# rCMAB = distance of each star from the CM of A and B
# vCMAB = speed of each star, w.r.t. the CM of A and B

# xvI_J = xv of star I with respect to the center of mass of star(s) J
# get output times
	xvB_A_AU		= AC.ReadAei(WhichDir, filenames[0], Bind[0], Bind[1])
	xvA_A_AU		= np.zeros_like(xvB_A_AU)
	if (mode == 'triple'):
		xvC_A_AU	= AC.ReadAei(WhichDir, filenames[1], Cind[0], Cind[1])
### Combine three stars into a 3D array
		xvA_AU=np.array([
                np.concatenate((xvA_A_AU,np.zeros((max(ntC-ntB,0.),6)))),
                np.concatenate((xvB_A_AU,np.zeros((max(ntC-ntB,0.),6)))),
                np.concatenate((xvC_A_AU,np.zeros((max(ntB-ntC,0.),6))))])
	else:
### Or two stars, if just a binary system
		xvA_AU=np.array([xvA_A_AU, xvB_A_AU])

### 1st dimension of array should be the number of stars
	assert(np.shape(xvA_AU)[0]) == nobjs
##################### Convert to mks units ##############################
	xvA = AC.AUtoMKS(xvA_AU)
### Get distances between stars (== rAij_AU*AU), and rel. velocity
	rAB = AC.Distance(xvA[0,:,:], xvA[1,:,:])
	vAB = AC.XVtoV(xvA[1,:,:])

######################## Binary system: #################################
# CMAB = Center of momentum frame of A+B
### Find the coordinates of the center of momentum
	xvCM_AB = AC.FindCM( m[0:2], xvA[0:2,0:ntB,:]) 
### Convert to center-of-momentum units
	xvAB  = AC.wrtCM(xvA[:,0:ntB,:], xvCM_AB)
### Get r, v in CM units
	rCMAB = np.array([ AC.XVtoR(xvAB[i,:,:]) for i in range(nobjs) ])
	vCMAB = np.array([ AC.XVtoV(xvAB[i,:,:]) for i in range(nobjs) ])
########################## Energies #####################################
# Get kinetic energy = (1/2)mv^2
	KCMAB = AC.Kinetic(m[0:2],vCMAB[0:2,:])
# Get potential energy = GMm/r
	UCMAB = np.array([ AC.Potential(m[0],m[1], rAB[0:ntB]),
					   AC.Potential(m[0],m[1], rAB[0:ntB]) ])
# Get total energy = K+U
	ECMAB = KCMAB+UCMAB

### This should be constant
	EtotCMAB = np.sum(ECMAB, 0)

### Get orbital parameters
	# grav. paramater for binary
	mu2 = G*(m[0]+m[1])
	# specific orbital energy
	epsB = AC.Eps(rAB[0:ntB], vAB[0:ntB], mu2)
	# semimajor axis
	aAB = AC.a(epsB, mu2)

	# specific angular momentum (r x v) of B wrt A
	hA     = AC.h( xvA[1,0:ntB,0:3],  xvA[1,0:ntB,3:6])
	hbarA  = AC.XVtoR( hA[:,:])
	# eccentricity
	eAB = AC.e(aAB, epsB, hbarA, mu2)
	# inclination
	iAB = AC.i(hA[:,2], hbarA)
	print(' aB = '+('% 10.4g'%(aAB[-1]/AU))
		 +', eB = '+('% 6.4g'%(eAB[-1]))
		 +', iB = '+('% 6.4g'%(iAB[-1])))

### Check for consistency
	dEpsB = epsB[-1]-epsB[0]

	if abs(dEpsB/epsB[0])<0.05:
		print(' dEpsB = '+ ('% 7.4g' % dEpsB)+' J, '+
			 ('% 7.4g' % float(100*dEpsB/epsB[0]))+'%'  )
	else:
		print(' dEpsB = '+('% 7.4g' % dEpsB)+' J, '+
			 ('% 7.4g' % float(100*dEpsB/epsB[0]))+'%'+
			  ' - large energy variation (AB)')

###################### Triple system: ########################################
### Only relevant if B and C both survived
	if (mode=='triple'):
#		nobjs = 3
		ind=min(ntB,ntC)
# get output times
#		LastxvB_A_AU = AC.ReadAei(WhichDir, filenames[0], ind-1, ind+4)
#		LastxvA_A_AU = np.zeros_like(LastxvB_A_AU)
#		LastxvC_A_AU = AC.ReadAei(WhichDir, filenames[1], ind-1, ind+4)

#		xvA_AU=np.array([
#                np.concatenate(( xvA_A_AU, np.zeros((max(ntC-ntB,0.),6)) )),
#                np.concatenate((xvB_A_AU,np.zeros((max(ntC-ntB,0.),6)))),
#                np.concatenate((xvC_A_AU,np.zeros((max(ntB-ntC,0.),6))))])

### If B is ejected before C, exend the AB CM, and set equal to A during that
		if ((ntB < ntC) & (DestB.strip(' ')=='ejected')):
			xvCM_AB = np.concatenate( (xvCM_AB[0:ntB,:], xvA[0,ntB:ntC,:]) )
		# Get xv wrt CM_AB through ntC
			xvAB  = AC.wrtCM(xvA[:,0:ntC,:], xvCM_AB)

### Combine three stars into a 3D array
#		xvA_Triple_AU = np.array([ xvA_A_AU[0:ind,:], 
#								   xvB_A_AU[0:ind,:], 
#								   xvC_A_AU[0:ind,:] ])
		xvA_Triple_AU = xvA_AU[:, 0:ntC, :]
### Convert to mks units
		xvA_Triple = AC.AUtoMKS(xvA_Triple_AU)
### suffix '2' => treating AB as one star at their CM, AB-C as binary
		xvA_Triple2 = np.array([ xvCM_AB[0:ntC,:], xvA[2,0:ntC,:] ])
### Find the coordinates of the center of momentum
		xvCM_ABC = AC.FindCM( m, xvA_Triple) 
		if ((ntB < ntC) & (DestB.strip(' ')=='ejected')):
			xvCM_ABC[0:ntB, :] = AC.FindCM( [m[0], 0., m[2]], xvA_Triple[:, 0:ntB, :]) 
#		xvCM_ABC2= AC.FindCM( [ m[0]+m[1], m[2] ], xvA_Triple2) 
### Convert to center-of-momentum units
		xvABC  = AC.wrtCM(xvA_Triple,  xvCM_ABC)
		xvABC2 = AC.wrtCM(xvA_Triple2, xvCM_ABC)
### Get r, v in CM units
		rCMABC = np.array([ AC.XVtoR( xvABC[i,:,:]) for i in range(nobjs) ])
		vCMABC = np.array([ AC.XVtoV( xvABC[i,:,:]) for i in range(nobjs) ])
		vCMABC2= np.array([ AC.XVtoV(xvABC2[i,:,:]) for i in range(2) ])
### Get distances between stars
		rAB_ABC = AC.Distance(xvABC[0,:,:], xvABC[1,:,:])
		rBC_ABC = AC.Distance(xvABC[1,:,:], xvABC[2,:,:])
		rAC_ABC = AC.Distance(xvABC[0,:,:], xvABC[2,:,:])
		rAB_C2  = AC.Distance(xvABC2[0,:,:], xvABC2[1,:,:])
########################## Energies #####################################
### Calculate kinetic energies, K
		KCMABC = AC.Kinetic(m, vCMABC)
		KABC2= AC.Kinetic( [m[0]+m[1], m[2]], vCMABC2)
		
### Calculate potential energies, U
		UAB_CMABC = AC.Potential(m[0],m[1], rAB_ABC)
		UBC_CMABC = AC.Potential(m[0],m[2], rAC_ABC)
		UAC_CMABC = AC.Potential(m[1],m[2], rBC_ABC)
		UCMABC = np.array([ UAB_CMABC+UAC_CMABC, 
                            UAB_CMABC+UBC_CMABC,
                            UBC_CMABC+UAC_CMABC ])
		UABC2    = np.array([ AC.Potential( m[0]+m[1], m[2], rAB_C2),
							  AC.Potential( m[0]+m[1], m[2], rAB_C2) ])
# Calculate total energy per object, E=K+U
		ECMABC = KCMABC+UCMABC
		ECMABC2= KABC2+UABC2
# Calculate total energy in the system
		EtotCMABC = np.sum(ECMABC, 0)
		EtotCMABC2= np.sum(ECMABC2, 0)
### Get orbital parameters
		# dist. from C to CM(AB)
		rC_AB = AC.XVtoR(xvAB[2,0:ntC,:])
		# vel. of C wrt CM(AB)
		vC_AB = AC.XVtoV(xvAB[2,0:ntC,:])
		# specific orbital energy
		epsC  = AC.Eps(rC_AB, vC_AB, mu)
		# semimajor axis
		aC = AC.a(epsC, mu)

		# specific angular momentum (r x v) of B wrt A
		hC    = AC.h(xvAB[2,0:(ntC),0:3], xvAB[2,0:(ntC),3:6])
		# magnitude of h
		hbarC = AC.XVtoR(hC[:,:])
		# eccentricity
		eC    = AC.e(aC, epsC, hbarC, mu)
		# inclination
		iC    = AC.i(hC[:,2], hbarC)
		print(' aC = '+('% 10.4g'%(aC[-1]/AU))
			 +', eC = '+('% 6.4g'%(eC[-1]))
			 +', iC = '+('% 6.4g'%(iC[-1])) )

### Check for consistency
		dEpsC = epsC[-1]-epsC[0]

		if abs(dEpsC/epsC[0])>0.02:
			print('dEpsC = {0} J, {1}% - large energy variation (ABC)'.format(
			     ('% 7.4g' % dEpsC), 
				 ('% 7.4g' % float(100*dEpsC/epsC[0])) ))
		else:
			print('dEpsC = {0} J'.format( ('% 7.4g' % dEpsC) ))

### Save dE for Prx sims
#		Esumnames = ['EBi','EBf','ECi','ECf','dEBe33','dECe33']
#		Esum = [   EtotCMAB[0],   EtotCMAB[-1], 
#				 EtotCMABC2[0], EtotCMABC2[-1], 
#					 dEAB/1e33,		dEABC/1e33]
#		if 'Prx' in WhichDir:
#			if '01' in WhichDir:
#				PrxSum=open('Proxlike/Plots/PrxSum.txt','w')
#				PrxSum.write(' '.join([('% 11s' % i) for i in Esumnames])+'\n')
#			else:
#				PrxSum=open('Proxlike/Plots/PrxSum.txt','a')
#			PrxSum.write(' '.join([('% 11.4g' % i) for i in Esum])+'\n')
#			PrxSum.close()

### Ejected objects should have pos. energy, stable ones should be neg.
	if ((DestC=='ejected') & (ECMABC[2,-1]<0.)):
		print('C: Energy weirdness!!!')
	if ((DestB=='ejected') & ( ECMAB[1,-1]<0.)):
		print('B: Energy weirdness!!!')

### Return all the data
	return 	  rAB,   vAB,  EtotCMAB,  ECMAB,  KCMAB,  UCMAB, dEpsB, epsB,\
			rC_AB, vC_AB, EtotCMABC, ECMABC, KCMABC, UCMABC, dEpsC, epsC,\
			aAB, eAB, iAB, aC, eC, iC, \
			UABC2, KABC2, ECMABC2, EtotCMABC2,	\
			t, ntB, ntC, ind, TimeB, DestB, TimeC, DestC

###############################################################################
def WriteAEI(WhichDir,ThisT,mode='triple'):
	'''Get the time-dependent data and write in a usable way to TimeData.txt'''
		
	print('WriteAEI      '+WhichDir)

### Needed modules
	import numpy as np
	import AlphaCenModule as AC
	from mks_constants import G, mSun, AU, day, m, mu

### Column width in output
	wn = [9]+2*[9,9,11,9,9,9]
	ws = [str(i) for i in wn]

### Get final orbit data from mercury's .aei files and analysis
	rB, vB,  EtotCMAB,  ECMAB,  KCMAB,  UCMAB, dEpsB, epsB,\
	rC, vC, EtotCMABC, ECMABC, KCMABC, UCMABC, dEpsC, epsC,\
	aB, eB, iB, aC, eC, iC,	\
	UABC2, KABC2, ECMABC2, EtotCMABC2, \
	t, ntB, ntC, ind, TimeB, DestB, TimeC, DestC = AC.GetFinalData(
													WhichDir, ThisT, mode)

### Make array of the binary and triple parameters over time
	data   = np.transpose(np.array([
				[(  '%9.3e' % i) for i in     t],
				[(  '%9.3f' % i) for i in rB/AU],
				[(  '%9.2f' % i) for i in    vB],
				[('% 11.4e' % i) for i in  epsB],
				[( '% 9.2f' % i) for i in aB/AU],
				[(  '%9.5f' % i) for i in    eB],
				[(  '%9.3f' % i) for i in    iB],
				[(  '%9.3f' % i) for i in rC/AU],
				[(  '%9.2e' % i) for i in    vC],
				[('% 11.4e' % i) for i in  epsC],
				[( '% 9.1f' % i) for i in aC/AU],
				[(  '%9.5f' % i) for i in    eC],
				[(  '%9.3f' % i) for i in    iC] ]))
	datalist = [' '.join(row)+'\n' for row in data]

	hdr = np.array(['t','rB', 'vB', 'epsB', 'aB', 'eB', 'iB',
						 'rC', 'vC', 'epsC', 'aC', 'eC', 'iC'])
	hdr = ' '.join([hdr[i].rjust(wn[i]) for i in range(len(hdr))])+'\n'

### Write data to file
	f = WhichDir+'/Out/AeiOutFiles/TimeData.txt'
	TimeFile=open(f,'w')
	TimeFile.write(hdr)
	for row in datalist:
		TimeFile.write(row)
	TimeFile.close()

###############################################################################
def Summary(WhichDir,ThisT,Tmax=1e9,WhichTime='1',machine='',
			wantsum=True,wantplot=False,mode='triple',cent='A'):
	'''A program to read the important bits and record in summary.txt'''
		
	print('	Summary      '+WhichDir+',                          WhichTime = '+
		  WhichTime)

### Needed modules
	import numpy as np
	import os
	import AlphaCenModule as AC
	from numpy import log10, sqrt, sin, pi
	from operator import add
	import subprocess
#	import rvtest as rv
	from mks_constants import G, mSun, AU, day, m, mu

	np.set_printoptions(precision=2)
### Stellar masses
	m = np.array(m)
	mA, mB, mC = m[0], m[1], m[2]

### Get initial parameters from 
	aeiIn=AC.InitParams(WhichDir)

### Get final orbits from element.out
#	B, C = AC.Elem(WhichDir)

### Get other final orbit data from .aei files and analysis
	rB, vB,  EtotCMAB,  ECMAB,  KCMAB,  UCMAB, dEpsB, epsB,\
	rC, vC, EtotCMABC, ECMABC, KCMABC, UCMABC, dEpsC, epsC,\
	aAB, eAB, iAB, aC, eC, iC,	\
	UABC2, KABC2, ECMABC2, EtotCMABC2, \
	t, ntB, ntC, ind, TimeB, DestB, TimeC, DestC = AC.GetFinalData(
													WhichDir, ThisT, mode)

##################################### Plot ################################
### Make plots of sim
	if ((machine != 'chloe') & (wantplot==True)):
		AC.MakePlots('binary',WhichDir, t, EtotCMAB, 
					 ECMAB, KCMAB, UCMAB, rB, suffix='_AB')
		if (mode=='triple'):
			AC.MakePlots('triple',WhichDir, t[0:ind], EtotCMABC2, 
						 ECMABC2, KABC2, UABC2, rC, suffix='_ABC')

#################################### Write ################################
### Arrange data nicely
	if ((wantsum==True) & (mode=='triple')):
		summaryfields=[aeiIn[0],aeiIn[1],aeiIn[2],aeiIn[3],aeiIn[4],aeiIn[5],
			    str(round(rB[-1]/AU,2)),
			    ('% 9.3g' % epsB[-1]),
			    str(round(aAB[-1]/AU,2)),
			    str(round(eAB[-1]   ,2)),
			    str(round(iAB[-1]   ,1)),
			    str(round(rC[ntC-1]/AU,1)),
			    ('% 9.3g' % epsC[-1]),
			    str(round( aC[-1]/AU,2)),
			    str(round( eC[-1]   ,2)),
			    str(round( iC[-1]   ,1)),
		        str(round(log10(float(TimeB)),5)), DestB,
			    str(round(log10(float(TimeC)),5)), DestC, 
				('% 9.3g' %dEpsB), ('% 9.3g' % dEpsC), '\n']
		headerfields=['aB','eB','iB','aC','eC','iC',
					 'rBf','EBf','aBf','eBf','iBf',
					 'rCf','ECf','aCf','eCf','iCf',
					 'logtB','destB','logtC','destC', 'dEAB','dEpsC', '\n']
		sumspaces=[
			len(aeiIn[0]),len(aeiIn[1]),len(aeiIn[2]),
				len(aeiIn[3]),len(aeiIn[4]),len(aeiIn[5]), 
			6,9,7,5,6, 9,9,9,6,6, 
			7,len(DestB),7,len(DestC), 9,9, len('\n')]
### Assemble summary and summary header rows with proper spacing

		summary=[summaryfields[i].rjust(sumspaces[i]) 
							for i in range(len(sumspaces))]
		header =[ headerfields[i].rjust(sumspaces[i]+1) 
							for i in range(len(sumspaces))]
		header[0]=header[0][2:]

#		header = '  aB    eB   iB     aC    eC    iC'+\
#			'    rBf       EBf'+\
#			'     aBf   eBf    iBf'+\
#			'       rCf       ECf'+\
#			'       aCf    eCf    iCf'+\
#			'   logtB    destB   logtC    destC\n'
#		summary = aeiIn+\
#			    [str(round(rAB[-1]/AU,2)).rjust(6)]			+\
#			    [('% 7.2g' % EtotCMAB[-1]).rjust(9)]		+\
#			    [str(round(aAB[-1]/AU,2)).rjust(7)]			+\
#			    [str(round(eAB[-1]   ,2)).rjust(5)]			+\
#			    [str(round(iAB[-1]   ,1)).rjust(6)]			+\
#			    [str(round(rCMAB[2,ntC-1]/AU,1)).rjust(9)]	+\
#			    [('% 7.2g' % EtotCMABC2[-1]).rjust(9)]		+\
#			    [str(round( aC[-1]/AU,2)).rjust(9)]			+\
#			    [str(round( eC[-1]   ,2)).rjust(6)]			+\
#			    [str(round( iC[-1]   ,1)).rjust(6)]			+\
#		        [str(round(log10(float(TimeB)),5)).rjust(7)]			+\
#				[DestB]										+\
#			    [str(round(log10(float(TimeC)),5)).rjust(7)]			+\
#				[DestC]+['\n']



### Determine if simulation is ending, and write data if so	
		AC.SummaryStatus(WhichDir, WhichTime, Tmax, ThisT, summary, header,
	                     rB/AU, EtotCMAB[:], rC/AU, EtotCMABC2[:],
	                     aAB[-1], eAB[-1], aC[-1], eC[-1], 
						 TimeB, TimeC, DestB, DestC)
	
############################################################################
def SummaryStatus(WhichDir, WhichTime, Tmax, ThisT, summary, summaryheader, 
                  rB, EB, rC, EC, 
                  aBf, eBf, aCf, eCf, 
				  TimeB, TimeC, DestB, DestC):
	'''Determine if simulation is ending, and write data if so	'''

### Needed modules
	import numpy as np
	import AlphaCenModule as AC
	from numpy import log10, sqrt, sin, pi
#	from operator import add
#	import subprocess
#	import rvtest as rv
	from mks_constants import AU

### Write to summary.out, but only if the simulation is ending:
### Determine status of each star in each survival criterion,
### i.e., whether an object is missing and the run can be ended
	isBmaxT = (float(TimeB)==Tmax)  # has either survived the max time?
	isCmaxT = (float(TimeC)==Tmax)	# or
	Bejectd = (DestB!='-'.rjust(8))	# has either been ejected/accreted?
	Cejectd = (DestC!='-'.rjust(8))	

### The stop conditions should be mutually exclusive
	if (isBmaxT==Bejectd==True) | (isCmaxT==Cejectd==True):
		print('Warning: possible stop condition conflict?')

### 'Little' stop = combination of these checks
	stop=(isBmaxT | isCmaxT | Bejectd | Cejectd)
	print('	    stop = '+
	str(isBmaxT)[0]+str(isCmaxT)[0]+str(Bejectd)[0]+str(Cejectd)[0]+
	',               Time = 1e'+str(int(log10(ThisT)))+' yrs')

### Write stop status to file for bash script to check
	StopFile=open(WhichDir+'/stopfile.txt','w')
	StopFile.write(str(stop)+'\n')
	StopFile.close()

### Prox's distance > this number (in AU) counts as 'prox-like'
	pcut=10000.
### If a Proxima-like C was created, stop the whole series of runs
	bigstop=False
	if (isBmaxT & isCmaxT):
		bigstop = ((aCf/AU)*(1+eCf) >= pcut) & (EC[-1] <= 0.)
### Weird circumstances that I want to stop and investigate:
	print('Testing for errors')
	tests=np.array([float(i) for i in summary[7:16]])
	if ( (float(EB[-1])<0.) & Bejectd ):
		if ( aBf*(1+eBf)>=1e5 ):
			print('B ejected due to extremely large orbit')
		else:
			bigstop = True
			print('**BIGSTOP FATE/ENERGY CONFLICT (B)**')			
	if ( (float(EC[-1])<0.) & Cejectd ):
		if ( aCf*(1+eCf)>=1e5 ):
			print('C ejected due to extremely large orbit')
		else:
			bigstop = True
			print('**BIGSTOP FATE/ENERGY CONFLICT (C)**')			
	if ((isBmaxT & isCmaxT) & (((float(EB[-1])>0.) & (not Bejectd)) | 
		 					   ((float(EC[-1])>0.) & (not Cejectd)) ) ):
		bigstop = True
		print('**BIGSTOP FATE/ENERGY CONFLICT**')			
	elif (np.isnan(float(EB[-1])) | np.isnan(float(EC[-1])) | 
	      np.isinf(float(EB[-1])) | np.isinf(float(EC[-1])) ):
		bigstop = True
		print('**BIGSTOP ENERGY ERROR**')
	elif ( any(np.isnan(tests))   | 
		   any(np.isinf(tests))   ):
		bigstop = True
		print('**BIGSTOP NONSENSICAL OUTPUTS**')
		print(tests)
	elif (Bejectd & Cejectd):
		bigstop = True
		print('**DOUBLE EJECTION -- TEST FOR CONSISTENCY**')
	else:
		print('  Error check passed')

	print('	 bigstop = '+str(bigstop).rjust(5)+',                       '+
		  '           WhichTime = '+str(WhichTime))

### debugging ----------
#	t1=float('-INF')
#	t2=float('-NAN')
#	print(t1,t2)
#	print(np.isnan(t1),np.isnan(t2))
#	print(np.isinf(t1),np.isinf(t2))
#	if (np.isnan(float(t1)) | np.isnan(float(t2)) | 
#	    np.isinf(float(t1)) | np.isinf(float(t2))):
#		print('true')
#	else:
#		print('false')
### debugging ----------

### Write big stop status to file for bash script to check
	BigStopFile=open(WhichDir+'/bigstopfile.txt','w')
	BigStopFile.write(str(bigstop)+'\n')
	BigStopFile.close()

### Write summary to file
	done=False
#	if "Prx" not in WhichDir:
	done=AC.WriteSummary(WhichDir,summary, summaryheader, stop, bigstop)

############################################################################
### Write summary to file
def WriteSummary(WhichDir, summary, summaryheader, stop, bigstop):

#	print('	WriteSummary '+WhichDir)

### Needed modules
	import os 
	import AlphaCenModule as AC

	sumpath=WhichDir+'/summary.out'

	if (stop==True | bigstop==True):
		SumFile=open(sumpath, 'a')
		if os.path.getsize(sumpath)==0:
			SumFile.write(''.join(summaryheader))
		strsum=" ".join(summary).strip()+'\n'
		SumFile.write(strsum)
		SumFile.close()

	return(True)

############################################################################
def MakePlots(version, WhichDir, t, Etot, E, K, U, r, suffix=''):
	'''Make plots of parameters from this run'''

	print('	MakePlots    '+WhichDir)

### Modules
	from mks_constants import G, mSun, AU, day, m, mu
	import numpy as np
	from numpy import sin, pi
	import matplotlib
	matplotlib.use('Agg', warn=False)
	import matplotlib.pyplot as plt

	nobj= E.shape[0]
	if (nobj != 2):
		print('Too many objects in E -- plot two objects')

# Assign colors
	if version=='binary':
		c=('b','y')
	elif version=='triple':
		c=('y','r')
	else:
		print('MakePlots version name invalid! Pick "binary" or "triple".')

# Plot energies vs. separation
	plt.plot(r/AU,   Etot, 'ks')
	for i in range(nobj):
		plt.plot(r/AU, E[i,:], c[i]+'s',
				 r/AU, U[i,:], c[i]+'o',
				 r/AU, K[i,:], c[i]+'^',
				)
#	plt.legend(('System total','A total','B total','A potential','A kinetic'),
#	           'lower right')
	plt.xlabel('Separation (AU)')
	plt.ylabel('Energy (J)')
	plt.title('Energies')
	plt.savefig(WhichDir+'/EvsR'+suffix+'.png')
	plt.clf()
	
# Plot energies over time
	plt.plot(t, Etot, 'k-')
	for i in range(nobj):
		plt.plot(t, E[i,:], c[i]+'-',
				 t, U[i,:], c[i]+'--',
				 t, K[i,:], c[i]+'-.',
	        )
	plt.xlabel('time (years)')
	plt.ylabel('Energy (J)')
	plt.title('Energies')
	plt.xscale('log')
	plt.savefig(WhichDir+'/EvsT'+suffix+'.png')
	plt.clf()
		
# Plot separation over time
	plt.plot(t, r/AU, 'k-')
	plt.xlabel('time (years)')
	plt.ylabel('Separation (AU)')
	plt.title('Distance over time')
	plt.xscale('log')
	plt.savefig(WhichDir+'/RvsT'+suffix+'.png')
	plt.clf()


# Plot C distance over time, if given
#	if not all( [value == 0. for value in rC] ):
#		plt.plot(t, rC/AU, 'r-')
##		for i in range(nobj):
##			plt.plot(t, E[i,:], c[i]+'-',
##					 t, U[i,:], c[i]+'--',
##					 t, K[i,:], c[i]+'-.',
##		        )
#		plt.xlabel('time (years)')
#		plt.ylabel('rC (AU)')
#		plt.title('C separation')
#		plt.xlim([9e8,10e8])
##		plt.xscale('log')
#		plt.savefig(WhichDir+'/C_RvsT'+suffix+'.png')
#		plt.clf()
## Plot C's energies over time
##		plt.plot(t, Etot, 'k-')
#		i=2
#		plt.plot(t, E[i,:], c[i]+'-',
#				 t, U[i,:], c[i]+'--',
#				 t, K[i,:], c[i]+'-.',
#	        )
#		plt.xlabel('time (years)')
#		plt.ylabel('Energy (J)')
#		plt.title('Energies')
##		plt.xscale('log')
#		plt.xlim([9e8,10e8])
#		plt.savefig(WhichDir+'/C_EvsT'+suffix+'.png')
#		plt.clf()
## Plot C's energies vs R
##		plt.plot(t, Etot, 'k-')
#		i=2
#		plt.plot(rC/AU, E[i,:], c[i]+'s',
#				 rC/AU, U[i,:], c[i]+'o',
#				 rC/AU, K[i,:], c[i]+'^',
#	        )
#		plt.xlabel('rC (AU)')
#		plt.ylabel('Energy (J)')
#		plt.title('Energies')
#		plt.savefig(WhichDir+'/C_EvsR'+suffix+'.png')
#		plt.clf()
## Plot C's energies vs R
##		plt.plot(t, Etot, 'k-')
#		plt.plot(rC/AU,   EtotC, 'ks')
#		for i in range(2):
#			plt.plot(rC/AU, EC[i,:], c[i+1]+'s',
#					 rC/AU, UC[:], c[i+1]+'o',
#					 rC/AU, KC[i,:], c[i+1]+'^',
#					)
#		plt.xlabel('rC (AU)')
#		plt.ylabel('Energy (J)')
#		plt.title('Energies')
#		plt.savefig(WhichDir+'/AB-C_EvsR'+suffix+'.png')
#		plt.clf()

############################################################################
def ReadInfo(WhichDir):
	'''A program to read the collision info from info.out'''

#	import numpy, os
#	import os 
	import AlphaCenModule as AC
	import numpy as np

 	InfoFile=open(WhichDir+'/Out/info.out','r')
	InfoLen=AC.FileLength(WhichDir+'/Out/info.out')
	AllInfo=InfoFile.readlines()
	InfoFile.close()

	start=where(AllInfo,"   Beginning the main integration.\n")
	ends=where(AllInfo,"   Integration complete.\n")
	nloops=len(ends)

	# Read the last loop only
#	if (nloops == 1):	
#		want=[AllInfo[j] for j in range(start[0]+2,ends[0]-1)]
#	else:
#		want=[AllInfo[j] for j in range(ends[nloops-2]+10,ends[nloops-1]-1)]
#	print(np.array(AllInfo))
#	print([AllInfo[j] for j in range(start[0]+2,ends[nloops-1]-1)])
#	print(start,ends)
#	print(nloops)

	# Get any line that isn't blank or boilerplate
	want = []
	for row in AllInfo[ (start[0]+2):(ends[nloops-1]-1) ]:
		if (row != '\n'):
			if not ( (row.split()[0] == 'Fractional' ) | 
					 (row.split()[0] == 'Integration') | 
					 (row.split()[0] == 'Continuing' ) ):
				want.append(row)
	
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
def SumAll(WhichDirs,cent,suffix=''):
	'''Read data from all directories and put in one file'''

#	import numpy, os, AlphaCenModule
	import os 
	import AlphaCenModule as AC
	
	print('	SumAll SumAll.out: '+", ".join(WhichDirs))

	Sum=[]
	Par=[]
	for j in range(len(WhichDirs)):
		DirSumPath=WhichDirs[j]+'/summary'+suffix+'.out'
		if os.path.isfile(DirSumPath):
### Read in summary.out from each directory, compile and write to SumAll.out
			DirSumFile=open(DirSumPath,'r')
			DirSumLen=AC.FileLength(DirSumPath)

			DirSum=DirSumFile.readlines()
			DirSum[0]='Dir '+DirSum[0]
			for k in range(1,len(DirSum)):
				DirSum[k]=str(j+1).rjust(3)+' '+DirSum[k]
			if (len(Sum) < 1):
				Sum=DirSum
			else:
				Sum=Sum+DirSum[1:]
			DirSumFile.close()
### Read in InParams.txt from each directory, compile into AllParams.txt
	 		DirParFile=open(WhichDirs[j]+'/InParams'+suffix+'.txt','r')
			DirParLen=AC.FileLength(WhichDirs[j]+'/InParams'+suffix+'.txt')
	
			DirPar=DirParFile.readlines()
			DirPar[0]='Dir '+DirPar[0]
			for k in range(1,len(DirPar)):
				DirPar[k]=str(j+1).rjust(3)+' '+DirPar[k]
			if (len(Par) < 1):
				Par=DirPar
			else:
				Par=Par+DirPar[1:]
			DirParFile.close()
	
### Write summary of summaries	
	SumAll=open('SumAll'+suffix+'.out','w')
	for j in range(len(Sum)):
		SumAll.write(Sum[j].rstrip()+'\n')
	SumAll.close()
### Write summary of initial parameters
	ParAll=open('AllParams'+suffix+'.txt','w')
	for j in range(len(Par)):
		ParAll.write(Par[j])
	ParAll.close()

############################################################################
#def InitSumFile(WhichDir):
#	'''Create an empty summary.out file if there isn't one'''
#
#### Needed modules
#	import os 
#
#	sumpath=WhichDir+'/summary.out'
#
#	SumFile=open(sumpath, 'a')
#	if os.path.getsize(sumpath)==0:
#		SumFile.write('    aB     eB     iB      aC     eC     iC'+\
#		'         rBf          EBf         rCf          ECf'+\
#		'            tB    destB            tC    destC\n')
#	SumFile.close()
###############################################################################





