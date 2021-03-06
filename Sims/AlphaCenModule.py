import AlphaCenModule as AC
import Merc
from mks_constants import G, mSun, mEarth, AU, day, deg2rad, rad2deg
import numpy as np
from numpy import pi, sin, cos, arccos, exp, log, log10, sqrt, linspace
from random import random, uniform
import os, re
import os.path
from operator import add
import subprocess
from operator import add
import socket
machine = socket.gethostname().split('.')[0]
if ('hammer' in machine):
	print("Warning: {0} doesn't have matplotlib!".format(machine))
else:
	from matplotlib import cm, colors
	import matplotlib.pyplot as plt

###############################################################################
def FileLength(fname):
	'''Function to count the number of lines in a file'''

	if os.path.getsize(fname)==0.:
		i = -2
	else:
#		with open(fname) as f:
		f = open(fname,'r')
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

###############################################################################
def WriteParam(WhichFile, stop, mA=1.105, step=10., user='yes', 
					output='default', dump='default'):
	'''Write param.in or param.dmp with specified parameters'''

	print('	WriteParam,          stop time = {0}'.format(stop))

	assert (user == 'yes') | (user == 'no')

### Read in the boilerplate text
	ParamBlankFile = open('paramblank.txt','r')
	ParamBlank = ParamBlankFile.readlines()
	ParamBlankFile.close()

	ParamText = np.array(ParamBlank)

### Assign values to any 'default' parameters
	NextT = np.ceil(np.log10(stop/365.25))
	if output == 'default':
		output =  365.25 * 10**(NextT - 3)
	if dump   == 'default':
		dump   = (365.25 * 10**NextT)/step/5

### Make lists of parameter info
	params = [stop, output, step, user, mA, int(dump)]
	inds   = [   7,      8,    9,   22, 28,   35]

### Insert parameters into the relevant lines
	for i,ind in enumerate(inds):
		ParamText[ind] = ParamText[ind].strip('\n')+str(params[i])+'\n'

### Write complete param text to file
	WriteParamFile = open(WhichFile,'w')
	for line in ParamText:
		WriteParamFile.write(line)
	WriteParamFile.close()

###############################################################################
def GetLastTime(WhichDir):
	'''Returns how far a simulation has gotten, based on info and aei.'''
	
### Modules

	AeiTime  = AC.GetAEILastTime(WhichDir)
	InfoTime = AC.GetInfoLastTime(WhichDir)
	
	LastTime = max([AeiTime,InfoTime])

	return LastTime

###############################################################################
def GetInfoLastTime(WhichDir):
	'''Returns how far a simulation has gotten, based on info.out.'''
	
	name,dest,time,PrevTime,complete = AC.ReadInfo(WhichDir)
	time = [float(i) for i in time]

### If Info was complete, final time = stoptime from param.dmp
	if (complete == True):
		ParamFile = open(WhichDir+'/Out/param.dmp','r')
		Param = ParamFile.readlines()
		ParamFile.close()
		for line in Param:
			if "stop time" in line:
				LastTime = float(line.split()[-1])/365.25
### If not complete, use max of last info starttime and collision times
	elif (complete == False):
		LastTime = max(time+[PrevTime])

	print('          Info time: {0}'.format(LastTime))
	return LastTime

###############################################################################
def GetAEILastTime(WhichDir):
	'''Returns end time of the last completed simulation step based on 
the .aei files. (Requires up-to-date element run for accuracy.)'''

	AeiDir = WhichDir+'/Out/AeiOutFiles/'
### Get list of files in AeiOutDir
	for root, dirs, files in os.walk(AeiDir):
		filelist=files
	filelist=np.array(filelist)
	validname = np.array([('.aei' in i) for i in filelist])
	filelist = filelist[validname]
### Get list of file sizes
	size=[]
	for f in filelist:
		size.append(os.path.getsize(AeiDir+f))
	size=np.array(size)
### Check which aei files are current (i.e. newer than the xv.out file)
	if os.path.isfile(WhichDir+'/Out/xv.out'):
		xvage = os.path.getmtime(WhichDir+'/Out/xv.out')
	elif os.path.isfile(WhichDir+'/Out/info.out'):
		print('No xv.out; using time from info.out to compare file ages')
		xvage = os.path.getmtime(WhichDir+'/Out/info.out')
	else:
		assert 0==1,'No xv or info.out files to compare ages!'
	age = []
	for f in filelist:
		age.append(os.path.getmtime(AeiDir+f)-xvage)
	age = np.array(age)
	isnew = (age >= 0.)
	assertfail='.aei files are older than xv.out file! Run element first'

### If none of the files are new, look at just the newest ones
	if not (sum(isnew)>=1):
		print("Warning! "+assertfail)
		mostrecent = max(age)
		isnew = np.array((age > (mostrecent-10.)))

### Get list of all files that are the maximum size (of the recent ones)
	isFullLength = np.array(size==max(size[isnew]))

### Want a file that is both recent and full length
	FullLengthFiles = filelist[isFullLength & isnew]

### Remove '.aei' from the end of the name of the first full-length
### file that is newer than the xv.out file
	FirstFileName = re.sub('\.aei$', '', FullLengthFiles[0])
### Get the maximum time 
	LastTime = AC.GetT(WhichDir,FirstFileName,-1,-0)[0]

	return(LastTime)

###########################################################################
def WriteObjInFile(WhichDir,names,filename,Header,FirstLines,xv,s,append='F'):
	'''Write big.in or small.in file'''

	if (append == 'F'):
		infile=open(WhichDir+'/In/'+filename+'.in','w')
	elif(append == 'T'):
		infile=open(WhichDir+'/In/'+filename+'.in','a')

### Header
	if (append == 'F'):
		for i in range(len(Header)):
			infile.write(Header[i])
### Data
	for i in range(len(names)):
		infile.write(FirstLines[i])
		infile.write("  {0: 24.20}  {1: 24.20}  {2: 24.20}\n".format(
												xv[i][0],xv[i][1],xv[i][2]))
		infile.write("  {0: 24.20}  {1: 24.20}  {2: 24.20}\n".format(
												xv[i][3],xv[i][4],xv[i][5]))
		infile.write("  {0: 24.20}  {1: 24.20}  {2: 24.20}\n".format(
												 s[i][0], s[i][1], s[i][2]))
	infile.close()

###############################################################################
def MakeBigRand(WhichDir,WhichTime, cent,
	aBmin, aBmax, eBmin, eBmax, iBmin, iBmax,
	aCmin, aCmax, eCmin, eCmax, iCmin, iCmax, 
	mA=1.105, mB=0.934, mC=0.123):
	'''Pick random parameters for the stars and make a new big.in'''

	print('	MakeBigRand  '+WhichDir+'/In/big.in,                          '+
		  str(WhichTime))

### Constants
#	mA, mB, mC = mA*mSun, mB*mSun, mC*mSun

### Pick random a, e, i and g, n, m for B and C
	aB = uniform(aBmin, aBmax)
	eB = uniform(eBmin, eBmax)
	iB = uniform(iBmin, iBmax)
	gB, nB, MB = uniform(0.0, 360.0), uniform(0.0, 360.0), uniform(0.0, 360.0)

	eC = uniform(eCmin, eCmax)
	iC = uniform(iCmin, iCmax)
	gC, nC, MC = uniform(0.0, 360.0), uniform(0.0, 360.0), uniform(0.0, 360.0)

	aCfactor = aB*(1+eB)/(1-eC)
	aC = uniform(aCmin*aCfactor, aCmax*aCfactor)

	aei = [[aB, eB, iB, gB, nB, MB], 
		   [aC, eC, iC, gC, nC, MC]]

### Read generic big.in file header; should move this into the main file	
	BigHeadFile=open('BigHeader.txt','r')
	BigHeader=BigHeadFile.readlines()
	BigHeadFile.close()

### First lines for each object
	if (cent == 'A'):
		BigFirstLines=([ 'AlCenB     m={0}  r=3.0\n'.format(mB)])
		names=['AlCenB']
	if (cent == 'B'):
		BigFirstLines=([ 'AlCenA     m={0}  r=3.0\n'.format(mA)])
		names=['AlCenA']
	BigFirstLines.append('PrxCen     m={0}  r=3.0\n'.format(mC))
	names.append('PrxCen')

### Spin
### No spin for all objects
	BigS=np.array([[0.0,  0.0,  0.0] for i in range(len(BigFirstLines))])

### Write big file
	AC.WriteObjInFile(
	WhichDir,names,'big',BigHeader,BigFirstLines,aei,BigS)

### Save initial parameters with full precision in InParams.txt
	InParams=open(WhichDir+'/InParams.txt','a')
	if os.path.getsize(WhichDir+'/InParams.txt')==0:
		InParams.write('                 aB                  eB'+\
		'                  iB                  aC                  eC'+\
		'                  iC                 gB                  nB'+\
		'                  MB                  gC                  nC'+\
		'                  MC\n')
	InParams.write(" ".join([repr(aB).rjust(19), repr(eB).rjust(19), 
													repr(iB).rjust(19), 
	repr(aC).rjust(19), repr(eC).rjust(19), repr(iC).rjust(19),
	repr(gB).rjust(19), repr(nB).rjust(19), repr(MB).rjust(19), 
	repr(gC).rjust(19), repr(nC).rjust(19), repr(MC).rjust(19),"\n"]))
	InParams.close()

###########################################################################
### Get orbital parameters of an object from big.in or small.in
def GetObjParams(filepath,obj):

### Read in *.in file
	f=open(filepath)
	infile=f.readlines()
	f.close()

### Find obj in infile
	for i,row in enumerate(infile):
		if (infile[i][0] != ')'):
			if (row.split()[0] == obj):
				objrow = i

### Extract orbital parameters (aei gnM or xyz uvw)
	x, y, z = infile[objrow+1].split()
	u, v, w = infile[objrow+2].split()
#	s1, s2, s3 = infile[objrow+3].split()

	return([float(x),float(y),float(z), float(u),float(v),float(w)])

###########################################################################
### Get orbital parameters of an object from big.in or small.in
def GetStellarMasses(WhichDir):

### Read in param.in and big.in files
	if   os.path.isfile(WhichDir+'/In/param.in'):
		par=open(WhichDir+'/In/param.in')
	elif os.path.isfile(WhichDir+'/Out/param.dmp'):
		par=open(WhichDir+'/Out/param.dmp')
	else:
		assert 0==1, 'No param.in or .dmp file found!'
	parfile=par.readlines()
	par.close()

	bigpath = WhichDir+'/In/big.in'
	if   os.path.isfile(bigpath):
		big=open(bigpath)
		nbiglines = FileLength(bigpath)
	elif os.path.isfile(WhichDir+'/Out/big.dmp'):
		big=open(WhichDir+'/Out/big.dmp')
		nbiglines = FileLength(WhichDir+'/Out/big.dmp')
	else:
		assert 0==1, 'No big.in or .dmp file found!'
	bigfile=big.readlines()
	big.close()
### Find central body mass
	for i,row in enumerate(parfile):
		if (row[0] != ')'):
			if 'central mass' in row:
				mstar = [float(row.split()[-1])]

### Identify header rows for each object in infile
	nbig = (nbiglines-6)/4
	objrows = [6+i*4 for i in range(nbig)]

### Extract masses
	for j,objrow in enumerate(objrows):
		massword = re.search('m *= *[0-9E+-.]*', bigfile[objrow]).group()
		mstar.append(float(massword.strip(' =m')))

	return(mstar)

###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   BUT WHAT FRAME ARE THEY IN???
def Merc_El2X(el, m):
	'''Convert orbital elements to cartesian for an ellipse (e < 1). 
Based on MCO_EL2X.FOR from J. Chambers' mercury6_2.for:
	gm = grav const * (central + secondary mass)
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of perihelion !!! )   in
	n = longitude of ascending node ) radians
	l = mean anomaly                )

	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(gm/a) )'''
# Big.in format:
# Asteroidal = Keplerian orbital elements, in an `asteroidal' format.
#              i.e.  a e I g n M, where
#               a = semi-major axis (in AU)
#               e = eccentricity
#               I = inclination (degrees)
#               g = argument of pericentre (degrees)
#               n = longitude of the ascending node (degrees)
#               M = mean anomaly (degrees)

### Extract needed parameters from input list
	a,e,i,g,n,M = [float(i) for i in el]
	gm = G*sum(m)

### Convert input degrees to radians
	i, g, n, M = deg2rad*i, deg2rad*g, deg2rad*n, deg2rad*M

### Pericenter
	q = (1.-e)*a

### Get p (longitude of pericenter) from g (argument of pericenter)
	p = g + n

### Rotation factors
	si, ci = sin(i), cos(i)
	sg, cg = sin(g), cos(g)
	sn, cn = sin(n), cos(n)
	z1 = cg * cn
	z2 = cg * sn
	z3 = sg * cn
	z4 = sg * sn
	d11 =  z1 - z4*ci
	d12 =  z2 + z3*ci
	d13 =       sg*si
	d21 = -z3 - z2*ci
	d22 = -z4 + z1*ci
	d23 =       cg*si

	assert (e < 1.0)	# haven't finished the other cases yet
### Calculate ellipse
	if (e < 1.0):
	### Ellipse
		romes = sqrt(1.0 - e*e)
		temp = AC.Merc_KeplerEllipse(e,M)
		se, ce = sin(temp), cos(temp)
		z1 = a * (ce - e)
		z2 = a * romes * se
		temp = sqrt(gm/a) / (1.0 - e*ce)
		z3 = -se * temp
		z4 = romes * ce * temp
#	elif (e == 1.0) then
	### Parabola
#		ce = orbel_zget(l)
#		z1 = q * (1.0 - ce*ce)
#		z2 = 2.0 * q * ce
#		z4 = sqrt(2.0*gm/q) / (1.0 + ce*ce)
#	 	z3 = -ce * z4
#	elif (e > 1.):
	### Hyperbola
#		romes = sqrt(e*e - 1.0)
#		temp = orbel_fhybrid(e,l)
#		se, ce = sin(temp), cos(temp)
#		z1 = a * (ce - e)
#		z2 = -a * romes * se
#		temp = sqrt(gm/abs(a)) / (e*ce - 1.0)
#		z3 = -se * temp
#		z4 = romes * ce * temp
	
###
	x = d11 * z1  +  d21 * z2
	y = d12 * z1  +  d22 * z2
	z = d13 * z1  +  d23 * z2
	u = d11 * z3  +  d21 * z4
	v = d12 * z3  +  d22 * z4
	w = d13 * z3  +  d23 * z4

	return(x,y,z,u,v,w)

###############################################################################
def Merc_KeplerEllipse(e,oldl):
	'''Solves Kepler's equation for eccentricities less than one.
 Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.

  e = eccentricity
  l = mean anomaly      (radians)
  u = eccentric anomaly (   "   )'''
#------------------------------------------------------------------------------

	twopi = 2.*pi
	piby2 = pi/2.

### Reduce mean anomaly to lie in the range 0 < l < pi
	if (oldl >= 0):
		l = oldl % twopi
	else:
		l = oldl % twopi + twopi
	sign = 1.0
	if (l > pi):
		l = twopi - l
		sign = -1.0

	ome = 1.0 - e
#-----------------------------
	if ((l >= .45) | (e < .55)):
### Regions A,B or C in Nijenhuis

### Rough starting value for eccentric anomaly
		if (l < ome):
			u1 = ome
		elif (l > (pi-1.0-e)):
			u1 = (l+e*pi)/(1.0+e)
		else:
			u1 = l + e

### Improved value using Halley's method
		flag = u1 > piby2
		if (flag):
			x = pi - u1
		else:
			x = u1
		x2 = x*x
		sn = x*(1.0 + x2*(-.16605 + x2*.00761) )
		dsn = 1.0 + x2*(-.49815 + x2*.03805)
		if (flag):
			dsn = -dsn
		f2 = e*sn
		f0 = u1 - f2 - l
		f1 = 1.0 - e*dsn
		u2 = u1 - f0/(f1 - .5*f0*f2/f1)

#---------------------
	else:
### Region D in Nijenhuis
### Rough starting value for eccentric anomaly
		z1 = 4.0*e + .50
		p = ome / z1
		q = .50 * l / z1
		p2 = p*p
		z2 = exp( log( sqrt( p2*p + q*q ) + q )/1.5 )
		u1 = 2.0*q / ( z2 + p + p2/z2 )

### Improved value using Newton's method
		z2 = u1*u1
		z3 = z2*z2
		u2 = u1 - .075*u1*z3 / (ome + z1*z2 + .375*z3)
		u2 = l + e*u2*( 3.0 - 4.0*u2*u2 )

### Accurate value using 3rd-order version of Newton's method
### N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!

### First get accurate values for u2 - sin(u2) and 1 - cos(u2)
	bigg = (u2 > piby2)
	if (bigg):
		z3 = pi - u2
	else:
		z3 = u2

	big = (z3 > (.5*piby2))
	if (big):
		x = piby2 - z3
	else:
		x = z3

	x2 = x*x
	ss = 1.0
	cc = 1.0

	ss = x*x2/6.*(1. - x2/20. *(1. - x2/42. *(1. - x2/72. *(1. -
         x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
	cc =   x2/2.*(1. - x2/12. *(1. - x2/30. *(1. - x2/56. *(1. -
         x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
         x2/306.))))))))

	if (big):
		z1 = cc + z3 - 1.0
		z2 = ss + z3 + 1.0 - piby2
	else:
		z1 = ss
		z2 = cc

	if (bigg):
		z1 = 2.0*u2 + z1 - pi
		z2 = 2.0 - z2

	f0 = l - u2*ome - e*z1
	f1 = ome + e*z2
	f2 = .5*e*(u2-z1)
	f3 = e/6.0*(1.0-z2)
	z1 = f0/f1
	z2 = f0/(f2*z1+f1)
	mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )

	return(mco_kep)

###############################################################################
### Add variations in pos/vel to a central/starting point



###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   for orbit with no inclination
def El2X(el, m):
	'''Convert orbital elements to cartesian for an ellipse (e < 1) with zero inclination
	mu = grav const * sum of masses
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of pericenter     )   in
	la = longitude of ascending node) radians
	f = true anomaly                )
	E = eccentric anomaly           )
	M = mean anomaly                )


	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(mu/a) )'''

### Extract needed parameters from input list
	a,e,i,g,la,f = [float(i) for i in el]
	mu = G*sum(m)

### Convert input degrees to radians
#	i, g, la = deg2rad*i, deg2rad*g, deg2rad*la
	f = deg2rad*f

### Calculate radius
	r = a*(1-e**2)/(1+e*cos(f))

### Position coords
	x = r*cos(f)
	y = r*sin(f)
	z = 0.

### Period
	T = ( (4*pi**2/mu) * a**3 )**0.5

### Mean motion
	n = 2*pi/T

### Velocity coords
	u = -    sin(f)  * n*a/(1-e**2)**0.5
	v =   (e+cos(f)) * n*a/(1-e**2)**0.5
	w = 0.

	return(x,y,z,u,v,w)
	
###########################################################################
def MakeSmallTestDisk(WhichDir,nmax=100,m='default',amin = 0.1,
	centobj='AlCenA', size='small'):
	'''Make a disk of small objects around A and B
 and write to small.in'''

	print('MakeSmallTestDisk {0}/In/small.in  n={1}'.format(
			WhichDir, nmax))

### Objects per disk (one disk each around A and B)
	nmax = int(nmax)
	num  = np.array(range(1,nmax+1))
	
### Other parameters
	if (m == 'default'):
		mtot = mEarth/1e6			# total disk mass
		m = mtot/nmax				# mass per disk particle
	r = 0.001						# Hill radii for small interactions
	d = 2.0							# small obj density, g/cm^3
	
### Get orbital parameters for central objects (in aei)
	mstars = np.array(AC.GetStellarMasses(WhichDir))*mSun

	aeiA, xvA = [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.]

	aeiB = AC.GetObjParams('{0}/In/big.in'.format(WhichDir), 'AlCenB')
	aeiB[0] = aeiB[0]*AU
	### Transform aei to xyz
	xvB  = AC.Merc_El2X(aeiB, [mstars[0],mstars[1]])
	xvB  = [xvB[0]/AU,     xvB[1]/AU,     xvB[2]/AU, 
			xvB[3]*day/AU, xvB[4]*day/AU, xvB[5]*day/AU]

### Small object parameters
	# List of small object names
	names    = ['M'+str(i) for i in num]

	# Determine disk boundaries (in AU)
	peri = aeiB[0]*(1-aeiB[1])/AU	# min distance btwn A and B
	amax = peri/2					# could experiment with this...
	# disk spacing = logorithmic
#	aspacing = amax - (amax-amin)*( log(nmax-n1+1)/log(nmax) )
	# or disk spacing = linear
	aspacing = (amin + (amax-amin)*( (num-1)/float(nmax-1) ))

	# aei parameters for each object, relative to central object
#	smallaei = np.array([[0. for j in range(6)] for i in n1])
	# Small object xv's
	SmlXV  =     np.array([[0. for j in range(6)] for i in num])
	if (size == 'big'):
		SmlAEI = np.array([[0. for j in range(6)] for i in num])
### Generate rocks
	for j in range(0,len(names)):
		
		# Determine parameters for central object
		if   (centobj=='AlCenA'):
			cent = xvA
			centmass = mstars[0]
		elif (centobj=='AlCenB'):
			cent = xvB
			centmass = mstars[1]
		else:
			print('what is going on? error with central object')

### Assign orbital parameters
		a        = aspacing[(num[j]-1)%nmax]*AU
		e,i, g,n = 0.0,0.0, 0.0,0.0,
		M        = 360*random()		# mean anomaly
#		E        =                  # eccentric anomaly
#		f        = 360*random()		# true anomaly
		aei = [a,e,i,g,n,M]

#		x,y,z, u,v,w = AC.El2X([a,e,i, g,n,f], [sum(mstars),m])
		x,y,z, u,v,w = AC.Merc_El2X(aei, [centmass,m])

		x,y,z, u,v,w = x/AU,y/AU,z/AU, u*day/AU,v*day/AU,w*day/AU
		xv = [x,y,z,u,v,w]

		### Coords = Jupiter/Saturn coords plus random variation
		obj = xv

		SmlXV[j]=[cent[i]+obj[i] for i in range(6)]

		if (size == 'big'):
			if   (j <  nmax):
				SmlAEI[j] = [aei[0]/AU,aei[1],aei[2],aei[3],aei[4],aei[5]]
			elif (j >= nmax):
				SmlAEI[j] = AC.Merc_X2El(SmallXV[j], [centmass,m])


### If needed to check velocities:
	vorb_actual = sqrt(SmlXV[:,3]**2+SmlXV[:,4]**2)
	vorb_expect = sqrt( G*(mstars[0]+m)/(aspacing*AU))*day/AU

#	print(np.transpose(np.array((vorb_actual,vorb_expect))))
#	matplotlib.use('Agg', warn=False)

#	plt.plot(aspacing, vorb_actual, 'b-',
#			 aspacing, vorb_expect, 'r--',
#			)
	#plt.legend(('Epsilon total','Potential','Kinetic'),
	#           'lower right')
	#plt.xlabel('time (years)')
	#plt.ylabel('Energy (J)')
	#plt.title('Energies')
	#plt.xscale('log')
#	plt.savefig('vorb_issue.png')
#	plt.clf()

### Read generic small.in file header; should move this into the main file	
	SmlHeadFile=open('SmlHeader.txt','r')
	SmlHeader=SmlHeadFile.readlines()
	SmlHeadFile.close()

### First lines for each object
	SmlParams="    m={0}  r={1} d={2}\n".format(m/mSun, r, d)
	SmlFirstLines=[]
	for i in names:
		SmlFirstLines.append(i+SmlParams)

### Spin
### No spin for all objects
	SmlS=np.array([[0.,0.,0.] for i in range(len(SmlFirstLines))])

### Write small file
	if (size=='small'):
		AC.WriteObjInFile(
		WhichDir,names,'small',SmlHeader,SmlFirstLines,SmlXV,SmlS)
	elif (size=='big'):
		AC.WriteObjInFile(
		WhichDir,names,'big',SmlHeader,SmlFirstLines,SmlAEI,SmlS,append='T')

###############################################################################
def InitParams(WhichDir):
	'''Read initial orbital parameters'''

### Digits
	iRnd =[1,3,1, 1,3,1]	# dec places to round starting a/e/i to, for B/C
	iJst = map(add, [4,2,3, 5,2,4], iRnd) # iRnd+iJst = # char per entry

### Get input parameters (aei for B and C)
	if os.path.isfile(WhichDir+'/InParams.txt') == False: 
		print('No InParams file!')
		aeiIn=['-','-','-','-','-','-']
	elif os.path.getsize(WhichDir+'/InParams.txt') == 0:
		print('Blank InParams file!')
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
def ReadAei(WhichDir, filename, index1=-1, index2=None):
	'''Read .aei file to get xyz and uvw (in AU and m/s respectively)'''

#	print('	--ReadAei      '+WhichDir+'/Out/AeiOutFiles/'+filename+', '+\
#		   str(index1)+':'+str(index2))


### Get last positions for surviving objects
	aeiFile=open(WhichDir+'/Out/AeiOutFiles/'+filename+'.aei','r')
	aei = aeiFile.readlines()
	if ((index2!=None) & (index1!=None)):
		xv = aei[index1:index2]
	elif ((index2==None) & (index1==None)):
		xv = aei[4:]
	elif (index2==None):
		xv = aei[index1:]
	elif (index1==None):
		xv = aei[4:index2]
### Arrange xyz uvw data into array
	xv = np.array([i.split()[6:12] for i in xv])
### Convert from strings to floats
	xv = xv.astype(np.float)

	return(xv)

###############################################################################
def GetT(WhichDir, filename, index1=-1, index2=0):
	'''Read .aei file to get time array. For whole file, index1=4'''

	print('	--GetT         '+WhichDir+'/Out/AeiOutFiles/'+filename+', '+\
		   str(index1)+':'+str(index2))

### Get last positions for surviving objects
	aeiFile=open(WhichDir+'/Out/AeiOutFiles/'+filename+'.aei','r')
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


	print('	--FindCM       m = '+str([('% 1.3g' % (i/mSun)) for i in m]))

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
#def Kinetic(m, v):
#	'''Calculate kinetic energy
#	where r is the distance between two stars, v is the relative velocity'''
#
#	K = [(1./2.)*m[i]*v[i,:]**2. for i in range(len(v[:,0]))]
#
#	return(np.array(K))
#
###############################################################################
#def Potential(m1, m2, r):
#	'''Calculate potential energy of one object due to another
#	where r is the distance between two stars, v is the relative velocity'''
#
#	U = -G*m1*m2/r
#	U = U/2.		# correction to avoid counting U twice
#
#	return(U)

###############################################################################
def Eps(r, v, mu):
	'''Calculate specific orbital energy
	where r is the distance between two stars, v is the relative velocity.
	Final units should be J/kg or m^2/s^2.'''

	eps = v**2./2.-mu/r

	return(eps)

###############################################################################
def kEps(v):
	'''Calculate specific orbital energy
	where r is the distance between two stars, v is the relative velocity.
	Final units should be J/kg or m^2/s^2.'''

	keps = v**2./2.

	return(keps)

###############################################################################
def uEps(r, mu):
	'''Calculate specific orbital energy
	where r is the distance between two stars, v is the relative velocity.
	Final units should be J/kg or m^2/s^2.'''

	ueps = -mu/r

	return(ueps)

###############################################################################
def mr(m):
	'''Calculate reduced mass'''

	m_1  = [1./i for i in m]
	mr_1 = sum(m_1)
	mr   = 1./mr_1

	return(mr)

###############################################################################
def h(r, v):
	'''Calculate angular momentum'''

	h = np.cross(r,v)

	return(np.array(h))

###############################################################################
def a(eps, mu):
	'''Calculate semimajor axis of an object's orbit'''

	a = -mu/(2.*eps)

	return(np.array(a))

###############################################################################
def e(eps, hbar, mu):
	'''Calculate eccentricity of an object's orbit'''

	e = (1.+2.*eps*hbar**2/mu**2)**0.5

	return(np.array(e))

###############################################################################
def i(hz, hbar):
	'''Calculate inclination of an object's orbit (in degrees)'''

	i = arccos(hz/hbar)*180./pi

	return(np.array(i))

###############################################################################
def GetFinalData(WhichDir,ThisT,mode, m, Tmax):
	'''Read in data, calculate derived values, and output for 
	writing/plotting/etc.'''

	print('	-GetFinalData  '+WhichDir)

	assert all(m > 10.)	# make sure units are kg, not mSun
	mu=G*sum(m)
	
	assert (mode=='triple') | (mode=='binary')

	if (mode == 'binary'):
		m=[m[0], m[1], 0.]
		filenames=['AlCenB']
	elif (mode == 'triple'):
		filenames=['AlCenB', 'PrxCen']

### Number of objects
	nobjs=len(filenames)+1
	print('            nobjs = {0}'.format(nobjs))
### Number of timesteps
	ntB=AC.FileLength(WhichDir+'/Out/AeiOutFiles/'+filenames[0]+'.aei')-4
	if (mode == 'triple'):
		ntC=AC.FileLength(WhichDir+'/Out/AeiOutFiles/'+filenames[1]+'.aei')-4
	else:
		ntC=0.
	if (ntB >= ntC):
		t = AC.GetT(WhichDir, filenames[0], 4, ntB+4)
	else:
		t = AC.GetT(WhichDir, filenames[1], 4, ntC+4)


	Bind = [4,ntB+4]
	Cind = [4,ntC+4]

### Fewest and most number of timesteps
	imin = min(ntB, ntC)
	imax = max(ntB, ntC)
### Default binary and triple system lengths
	iBin = imax
	iTri = imin

#------------------------------------------------------------------------------
########## Get object's fate and collision/ejection time from info.out #######
	name,dest,time,LastTime,complete = AC.ReadInfo(WhichDir)
	DestB,TimeB = '-'.rjust(8),str(ThisT).rjust(13)
	DestC,TimeC = '-'.rjust(8),str(ThisT).rjust(13)
	for j in range(len(name)):
		if (name[j] == 'AlCenB'):
			DestB = dest[j].rjust(8)
			TimeB = time[j].rjust(13)
		if (name[j] == 'PrxCen'):
			DestC = dest[j].rjust(8)
			TimeC = time[j].rjust(13)

### Get flags for various possible fatess
	isBmaxT, isCmaxT, Bejectd, Cejectd, \
	Bcollsn, Ccollsn, Bremovd, Cremovd = AC.ParseDestinations(DestB, DestC, \
														TimeB, TimeC, Tmax)

### If collision, have to add masses
	if (Bcollsn):
		if (DestB.strip() == 'Center'):
			m[0:2] = [m[0]+m[1], 0.]
		if (DestB.strip() == 'PrxCen'):
			m[1:3] = [0., m[0]+m[1]]
	if (Ccollsn):
		if (DestC.strip() == 'Center'):
			m[0:3] = [m[0]+m[1], m[2], 0.]
		if (DestC.strip() == 'AlCenB'):
			m[1:3] = [m[0]+m[1], 0.]
#		mu should be unchanged; mu2 is calculated later

#------------------------------------------------------------------------------
############# Read in original x and v values ###########################
#------------------------------------------------------------------------------
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
		nanCB, nanBC =  np.empty((max(ntC-ntB,0.),6)), \
						np.empty((max(ntB-ntC,0.),6))
		nanCB[:], nanBC[:] = np.NAN, np.NAN
		xvC_A_AU	= AC.ReadAei(WhichDir, filenames[1], Cind[0], Cind[1])
### Combine three stars into a 3D array
		xvA_AU=np.array([
                np.concatenate( (xvA_A_AU, nanCB) ),
                np.concatenate( (xvB_A_AU, nanCB) ),
                np.concatenate( (xvC_A_AU, nanBC) )])
	else:
### Or two stars, if just a binary system
		xvA_AU=np.array([xvA_A_AU, xvB_A_AU])
	if ((any(np.isnan(xvB_A_AU[-1,:]))) | (any(np.isinf(xvB_A_AU[-1,:])))):
		MercNanError = True
	elif ( mode=='triple'):
		if (any(np.isnan(xvC_A_AU[-1,:])) | any(np.isinf(xvC_A_AU[-1,:]))):
			MercNanError = True
#		if ((not any(np.isnan(np.concatenate((xvB_A_AU[-2,:],xvC_A_AU[-2,:]))))) |
#		    (not any(np.isinf(np.concatenate((xvB_A_AU[-2,:],xvC_A_AU[-2,:])))))):
#			print('NaNs in last timestep only; backing up by one')
#			iBin, iTri = iBin-1, iTri-1
		# placeholder so code runs even though I haven't fixed the above clause
		else:
			MercNanError = False
	else:
		MercNanError = False

### 1st dimension of array should be the number of stars
	assert(np.shape(xvA_AU)[0]) == nobjs
### Convert to mks units
	xvA = AC.AUtoMKS(xvA_AU)

### Get distances between stars (== rAij_AU*AU), and rel. velocity
	rAB = AC.Distance(xvA[0,:,:], xvA[1,:,:])
	if (mode=='triple'):
		rBC = AC.Distance(xvA[1,:,:], xvA[2,:,:])
		rAC = AC.Distance(xvA[0,:,:], xvA[2,:,:])


### Determine system classification based on last triple timestep
	if (mode=='triple'):
		distances = np.array([rAB[iTri-1], rBC[iTri-1], rAC[iTri-1]])
		if (((DestB.strip() != 'Center') & (DestB.strip() != 'PrxCen')) & 
			((DestC.strip() != 'Center') & (DestC.strip() != 'AlCenB'))):
			if   ( (min(distances) == rAB[iTri-1]) | (ntB > ntC)):
				sysclass = 'AB'
				abc = [0,1,2]
				if (ntB  < ntC):
					iBin = imin
			elif ((min(distances) == rAC[iTri-1]) | (ntB < ntC)):
				sysclass = 'AC'
				abc = [0,2,1]
				if (ntB  > ntC):
					iBin = imin
			elif (min(distances) == rBC[iTri-1]):
				sysclass = 'BC'
				abc = [1,2,0]
				if (ntB != ntC):
					iBin = imin
			else:
				print('no class assigned...???')
		else:
			# If any objects collide, assume AB binary?
			sysclass = 'AB'
			abc = [0,1,2]
	else:
		# If any objects collide, assume AB binary?
		sysclass = 'AB'
		abc = [0,1,float('NaN')]
		

### Print basic info about simulation
	print('       {4}, {0} ({2:6.4g} y), {1} ({3:6.4g} y)'.format(
			       DestB,       DestC, 
			float(TimeB),float(TimeC), sysclass))

# Get binary data
	aB, eB, iB, epsB, xvCM_AB, kB, uB, rB, vB = AC.Binary(
			[ m[abc[0]],m[abc[1]] ], 
			np.array([xvA[abc[0],:,:], xvA[abc[1],:,:]]), 
			iBin)
	xvAB = AC.wrtCM(np.array([xvA[abc[0],:,:], xvA[abc[1],:,:]]), xvCM_AB)
# Get triple system data
	if (mode=='triple'):
		aC, eC, iC, epsC, kC, uC, rC_AB, vC_AB = AC.Triple(
			[ m[abc[0]],m[abc[1]],m[abc[2]] ], 
			np.array([xvA[abc[0],:,:], xvA[abc[1],:,:], xvA[abc[2],:,:]]), 
			iTri, DestB, xvA_AU, xvCM_AB, 
		AC.wrtCM(np.array([xvA[abc[0],:,:], xvA[abc[1],:,:],xvA[abc[2],:,:]]), xvCM_AB))
#------------------------------------------------------------------------------
### Print results
	print('	   aB = {0:10.4g}, eB = {1:6.4g}, iB = {2:6.4g}'.format(
			   aB[-1]/AU,  eB[-1],    iB[-1]         ) )
	if (mode=='triple'):
		print('	   aC = {0:10.4g}, eC = {1:6.4g}, iC = {2:6.4g}'.format(
				   aC[-1]/AU,  eC[-1],    iC[-1]         ) )

### Check for consistency
### Ejected objects should have pos. energy, stable ones should be neg.
	dEpsB = epsB[-1]-epsB[0]
	if ((DestB=='ejected') & (epsB[-1]<0.)):
		print('B: Energy weirdness!!!')
	if (mode=='triple'):
		dEpsC = epsC[-1]-epsC[0]
		if ((DestC=='ejected') & (epsC[-1]<0.)):
			print('C: Energy weirdness!!!')

#------------------------------------------------------------------------------
	if (mode=='binary'):
		rC_AB,vC_AB, dEpsC, epsC, kC, uC,aC, eC, iC = [float('NaN')]*9
#		print(TimeC, DestC)
### Return all the data
	return	m, iBin, iTri, imin, imax, MercNanError,\
	   rB,    vB, dEpsB, epsB, kB, uB, \
	rC_AB, vC_AB, dEpsC, epsC, kC, uC, \
	aB, eB, iB, aC, eC, iC,	\
	t, TimeB, DestB, TimeC, DestC

#------------------------------------------------------------------------------
######################## Binary system: #################################
#------------------------------------------------------------------------------
def Binary(m, xv, ind):
	'''Calculate orbital parameters for the binary of the system, stars 1 & 2.
 Usually the binary will be stars A and B, but substitute the values of A & C 
or B & C if they are actually closer. They'll still be referred to as A & B
 here.'''

### Relative distance/velocity of binary stars
	r2 = AC.Distance(xv[0,0:ind,:], xv[1,0:ind,:])
	v2 =    AC.XVtoV(xv[1,0:ind,:]- xv[0,0:ind,:])


# CMAB = Center of momentum frame of A+B
### Find the coordinates of the center of momentum
	xvCM_2 = AC.FindCM( m[0:2], xv[0:2,0:ind,:]) 
### Convert to center-of-momentum units
	xv2  = AC.wrtCM(xv[:,0:ind,:], xvCM_2)
### Get r, v in CM units
	rCM2 = np.array([ AC.XVtoR(xv2[i,:,:]) for i in range(2) ])
	vCM2 = np.array([ AC.XVtoV(xv2[i,:,:]) for i in range(2) ])
#---------------------------- Energies ----------------------------------------
### Get orbital parameters
	# grav. paramater for binary
	mu2  = G*sum(m)
	# specific orbital energy
	eps2 = AC.Eps(r2[0:ind], v2[0:ind], mu2)
	k2   = AC.kEps(v2[0:ind])
	u2   = AC.uEps(r2[0:ind], mu2)
	# semimajor axis
	a2 = AC.a(eps2, mu2)

	# specific angular momentum (r x v) of B wrt A
	h2     = AC.h( xv[1,0:ind,0:3]-xv[0,0:ind,0:3],  
				   xv[1,0:ind,3:6]-xv[0,0:ind,3:6])
	hbar2  = AC.XVtoR( h2[:,:])
	# eccentricity
	e2 = AC.e(eps2, hbar2, mu2)
	# inclination
	i2 = AC.i(h2[:,2], hbar2)

	return(a2, e2, i2, eps2, xvCM_2, k2, u2, r2, v2)

#------------------------------------------------------------------------------
####################### Triple system: ########################################
#------------------------------------------------------------------------------
def Triple(m, xv, ind, DestB, xvA_AU, xvCM_AB, xvAB):
### 
	mu = G*sum(m)
### Combine xv data into 3D array in binary-triple order (usually ABC):
#	xv = np.array([xv1,xv2,xv3])

### If B is ejected before C, extend the AB CM, and set equal to A during that
#	if ((ntB < ntC) & (DestB.strip(' ')=='ejected')):
#		xvCM_AB = np.concatenate( (xvCM_AB[0:ntB,:], xv[0,ntB:ntC,:]) )
	# Get xv wrt CM_AB through ntC
#		xvAB  = AC.wrtCM(xv[:,0:ntC,:], xvCM_AB)

#----------------------------------- Energies ---------------------------------
### Get orbital parameters
	# dist. from C to CM(AB)
	rC_AB = AC.XVtoR(xvAB[2,0:ind,:])
	# vel. of C wrt CM(AB)
	vC_AB = AC.XVtoV(xvAB[2,0:ind,:])
### Specific orbital energy
	epsC  = AC.Eps(rC_AB, vC_AB, mu)
	kC    = AC.kEps(vC_AB)
	uC    = AC.uEps(rC_AB, mu)
### Orbital parameters
	# semimajor axis
	aC = AC.a(epsC, mu)
	# specific angular momentum (r x v) of B wrt A
	hC    = AC.h(xvAB[2,0:ind,0:3], xvAB[2,0:ind,3:6])
	# magnitude of h
	hbarC = AC.XVtoR(hC[:,:])
	# eccentricity
	eC    = AC.e(epsC, hbarC, mu)
	# inclination
	iC    = AC.i(hC[:,2], hbarC)
	return(aC, eC, iC, epsC, kC, uC, rC_AB, vC_AB)

###############################################################################
def WriteTimeData(WhichDir,ThisT='dflt',m='dflt',mode='dflt',Tmax='dflt'):
	'''Get the time-dependent data and write in a usable way to TimeData.txt'''
		
	print('        --WriteTimeData      '+WhichDir)

### Calculate default values from data in directory
	if m == 'dflt':
		m = AC.GetStellarMasses(WhichDir)
	assert ( (len(m) == 2) or (len(m)==3) ), 'Wrong # of objs -- designed for 2 or 3!'
	m  = np.array([i*mSun for i in m])
	mu = G*sum(m)

	if Tmax == 'dflt':
		Tmax = AC.GetLastTime(WhichDir)
	if ThisT == 'dflt':
		ThisT = Tmax
	
	if mode == 'dflt':
		if '-2' in WhichDir:
			mode = 'binary'
		else:
			mode = 'triple'
### Column width in output
	wn = [9]+2*[9,9,11,9,9,9,9]+[9]
	ws = [str(i) for i in wn]

### Get final orbit data from mercury's .aei files and analysis
	m, iBin, iTri, imin, imax, MercNanError,\
	rB, vB, dEpsB, epsB, kB, uB, \
	rC, vC, dEpsC, epsC, kC, uC, \
	aB, eB, iB, aC, eC, iC,	\
	t, TimeB, DestB, TimeC, DestC = AC.GetFinalData(
											WhichDir, ThisT, mode, m, Tmax)

	if (mode == 'triple') and (len(aC) < len(aB)):
		diff = len(rB)-len(rC)
		rC   = np.concatenate((  rC,[float('NaN') for i in range(diff)]))
		vC   = np.concatenate((  vC,[float('NaN') for i in range(diff)]))
		epsC = np.concatenate((epsC,[float('NaN') for i in range(diff)]))
		aC   = np.concatenate((  aC,[float('NaN') for i in range(diff)]))
		eC   = np.concatenate((  eC,[float('NaN') for i in range(diff)]))
		iC   = np.concatenate((  iC,[float('NaN') for i in range(diff)]))

	if mode == 'binary':
		testnans = [np.isnan(rC), np.isnan(vC), np.isnan(epsC), np.isnan(aC), np.isnan(eC), np.isnan(iC)]
		assert all(testnans), 'Binary mode but C values not NaN: {}'.format(testnans)
		rC   = np.array([float('NaN') for i in range(len(rB))])
		vC   = np.array([float('NaN') for i in range(len(rB))])
		epsC = np.array([float('NaN') for i in range(len(rB))])
		aC   = np.array([float('NaN') for i in range(len(rB))])
		eC   = np.array([float('NaN') for i in range(len(rB))])
		iC   = np.array([float('NaN') for i in range(len(rB))])

	pB = aB*(1-eB)
	pC = aC*(1-eC)

### If disk was simulated, find truncation radius
	if 'Disk' in WhichDir:
		dummy,dummy,dummy,rtr = AC.CalcDisk(WhichDir)
		assert len(rtr)>=len(rB),'Disk aei files not full length -- outdated?'
		if len(rtr)>len(rB):
			print('Warning: stars live shorter than disk particles; truncating rtr')
			rtr = rtr[:len(rB)]
	else:
		rtr = np.array([float('NaN') for i in range(len(rB))])

### Make array of the binary and triple parameters over time
	data   = np.transpose(np.array([
				[(  '%9.3e' % i) for i in     t],
				[(  '%9.3f' % i) for i in rB/AU],
				[(  '%9.2f' % i) for i in    vB],
				[('% 11.4e' % i) for i in  epsB],
				[( '% 9.2f' % i) for i in aB/AU],
				[(  '%9.5f' % i) for i in    eB],
				[(  '%9.3f' % i) for i in    iB],
				[( '% 9.2f' % i) for i in pB/AU],
				[(  '%9.3f' % i) for i in rC/AU],
				[(  '%9.2e' % i) for i in    vC],
				[('% 11.4e' % i) for i in  epsC],
				[( '% 9.1f' % i) for i in aC/AU],
				[(  '%9.5f' % i) for i in    eC],
				[(  '%9.3f' % i) for i in    iC],
				[( '% 9.2f' % i) for i in pC/AU],
				[( '% 9.2f' % i) for i in   rtr] ]))
	datalist = [' '.join(row)+'\n' for row in data]

	hdr = np.array(['t','rB', 'vB', 'epsB', 'aB', 'eB', 'iB','pB',
						 'rC', 'vC', 'epsC', 'aC', 'eC', 'iC','pC','rtr'])
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
			wantsum=True,wantplot=False,mode='triple',cent='A',
			mA=1.105, mB=0.934, mC=0.123):
	'''A program to read the important bits and record in summary.txt'''
		
	print('	Summary        '+WhichDir+',                          WhichTime = '+
		  WhichTime)

	np.set_printoptions(precision=2)
### Stellar masses
	mA, mB, mC = mA*mSun, mB*mSun, mC*mSun
	m  = np.array([mA,mB,mC])
	M  = sum(m)
	mu = G*M

### Get initial parameters from 
	aeiIn=AC.InitParams(WhichDir)

### Get final orbits from element.out
#	B, C = AC.Elem(WhichDir)

### Get other final orbit data from .aei files and analysis
	m, iBin, iTri, imin, imax, MercNanError,\
	rB, vB, dEpsB, epsB, kB, uB, \
	rC, vC, dEpsC, epsC, kC, uC, \
	aB, eB, iB, aC, eC, iC,	\
	t, TimeB, DestB, TimeC, DestC = AC.GetFinalData(
											WhichDir, ThisT, mode, m, Tmax)

##################################### Plot ################################
### Make plots of sim
	if (wantplot==True):
		if (machine != 'chloe'):
			AC.MakePlots(	 'binary',WhichDir, t[0:iBin], epsB, kB, uB,
						 		rB, m, suffix='_AB')
			if (mode=='triple'):
				AC.MakePlots('triple',WhichDir, t[0:iTri], epsC, kC, uC,
								rC, m, suffix='_ABC')
# Plot energies over time
			plt.plot(t[0:iBin], epsB*(m[0]+m[1]),      'r-')
			plt.plot(t[0:iTri], epsC*(m[0]+m[1]+m[2]), 'b-')
			plt.xlabel('time (years)')
			plt.ylabel('Epsilon')
			plt.title('Energy')
			plt.xscale('log')
			plt.legend(('Binary epsilon','Triple epsilon'),
	           'lower right')
			plt.savefig(WhichDir+'/EpsvsT.png')
			plt.clf()

			plt.plot(t[0:iBin], aB/AU, 'r-')
			plt.plot(t[0:iTri], aC/AU, 'b-')
			plt.xlabel('time (years)')
			plt.ylabel('a')
			plt.title('Semimajor Axis')
			plt.xscale('log')
			plt.yscale('log')
			plt.legend(('Binary','Triple'),
	           'upper left')
			plt.savefig(WhichDir+'/AvsT.png')
			plt.clf()

			plt.plot(t[0:iBin], eB, 'r-')
			plt.plot(t[0:iTri], eC, 'b-')
			plt.xlabel('time (years)')
			plt.ylabel('e')
			plt.title('Eccentricity')
			plt.xscale('log')
			plt.legend(('Binary','Triple'),
	           'upper left')
			plt.savefig(WhichDir+'/EvsT.png')
			plt.clf()

			plt.plot(aB/AU,eB, 'r-')
			plt.plot(aC/AU,eC, 'b-')
			plt.xlabel('a')
			plt.ylabel('e')
			plt.title('A vx E')
			plt.xscale('log')
			plt.legend(('Binary','Triple'),
	           'lower right')
			plt.savefig(WhichDir+'/AvsE.png')
			plt.clf()
##################################### Write ################################
### Arrange data nicely
	if ((mode=='triple')):
		summaryfields=[aeiIn[0],aeiIn[1],aeiIn[2],aeiIn[3],aeiIn[4],aeiIn[5],
			    str(round(rB[-1]/AU,2)),
			    ('% 9.3g' % epsB[-1]),
			    str(round(aB[-1]/AU,2)),
			    str(round(eB[-1]   ,2)),
			    str(round(iB[-1]   ,1)),
			    str(round(rC[iTri-1]/AU,1)),
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

### Determine if simulation is ending, and write data if so	
		AC.SummaryStatus(WhichDir, WhichTime, Tmax, ThisT, summary, header,
	                     wantsum, MercNanError, rB/AU, epsB, rC/AU, epsC,
	                     aB[-1], eB[-1], aC[-1], eC[-1], 
						 TimeB, TimeC, DestB, DestC, iBin, imax)
	
############################################################################
def ParseDestinations(DestB, DestC, TimeB, TimeC, Tmax):
	'''Determine what hit what, etc.'''

### Write to summary.out, but only if the simulation is ending:
### Determine status of each star in each survival criterion,
### i.e., whether an object is missing and the run can be ended
	isBmaxT = (float(TimeB)==Tmax)  # has either survived the max time?
	isCmaxT = (float(TimeC)==Tmax)	# or
	Bejectd = (DestB.strip()=='ejected')	# has either been ejected/accreted?
	Cejectd = (DestC.strip()=='ejected')
	Bcollsn = ((DestB.strip()!='-') & (DestB.strip()!='ejected'))
	Ccollsn = ((DestC.strip()!='-') & (DestC.strip()!='ejected'))	
	Bremovd = Bejectd | Bcollsn
	Cremovd = Cejectd | Ccollsn

	return(	isBmaxT, isCmaxT, \
			Bejectd, Cejectd, \
			Bcollsn, Ccollsn, \
			Bremovd, Cremovd)

############################################################################
def SummaryStatus(WhichDir, WhichTime, Tmax, ThisT, summary, summaryheader, 
                  wantsum, MercNanError, rB, epsB, rC, epsC, 
                  aBf, eBf, aCf, eCf, 
				  TimeB, TimeC, DestB, DestC, iBin, imax):
	'''Determine if simulation is ending, and write data if so	'''


### Get flags for various possible fatess
	isBmaxT, isCmaxT, Bejectd, Cejectd, \
	Bcollsn, Ccollsn, Bremovd, Cremovd = AC.ParseDestinations(DestB, DestC,
															   TimeB, TimeC,
																		 Tmax)

### The stop conditions should be mutually exclusive
	if (isBmaxT==Bremovd==True) | (isCmaxT==Cremovd==True):
		print('Warning: possible stop condition conflict?')

### 'Little' stop = combination of these checks
	stop=(isBmaxT | isCmaxT | (Bremovd) | (Cremovd))
	print('	    stop = '+
	str(isBmaxT)[0]+str(isCmaxT)[0]+str(Bremovd)[0]+str(Cremovd)[0]+
	',               Time = 1e'+str(int(log10(ThisT)))+' yrs')

### Write stop status to file for bash script to check
	if (wantsum == True):
		StopFile=open(WhichDir+'/stopfile.txt','w')
		StopFile.write(str(stop)+'\n')
		StopFile.close()

### Prox's distance > this number (in AU) counts as 'prox-like'
	pcut=10000.
### If a Proxima-like C was created, stop the whole series of runs
	bigstop=False
	err=False
	if (isBmaxT & isCmaxT):
		bigstop = ((aCf/AU)*(1+eCf) >= pcut) & (epsC[-1] <= 0.)
### Weird circumstances that I want to stop and investigate:
	print('	Testing for errors')
	tests=np.array([float(i) for i in summary[7:16]])
	if ( (float(epsB[-1])<0.) & Bejectd ):
		if ( aBf*(1+eBf)/AU >= 1e5 ):
			print('B ejected due to extremely large orbit')
		else:
			bigstop = True
			err=True
			print('**B EJECTED WHEN "STABLE"!**')			
	if ( (float(epsC[-1])<0.) & Cejectd ):
		if ( aCf*(1+eCf)/AU >= 1e5 ):
			err=True
			print('C ejected due to extremely large orbit')
		else:
			bigstop = True
			err=True
			print('**C EJECTED WHEN "STABLE"!**')			
	if ((isBmaxT & isCmaxT) & (((float(epsB[-1])>0.) & (not Bremovd)) | 
		 					   ((float(epsC[-1])>0.) & (not Cremovd)) ) ):
		bigstop = True
		err=True
		print('**STAR REMAINED WHEN "UNSTABLE"!**')			
	elif (np.isnan(float(epsB[-1])) | np.isnan(float(epsC[-1])) | 
	      np.isinf(float(epsB[-1])) | np.isinf(float(epsC[-1])) ):
		if MercNanError:
			err=True
			print('	 Mercury write error -- continuing')
		else:
			bigstop = True
			err=True
			print('**BIGSTOP ENERGY ERROR**')
	elif ( any(np.isnan(tests))   | 
		   any(np.isinf(tests))   ):
		bigstop = True
		err=True
		print('**BIGSTOP NONSENSICAL OUTPUTS**')
		print(tests)
	elif (iBin < imax):
		if (aBf < 0.):
			print('	 Binary system breakup -- continuing')
		else:
			bigstop = True
			print('**BINARY BREAKUP WHEN "STABLE"**')
	else:
		print('	  Error check passed')
		if ((bigstop == True) & (err == False)):
			print('   => *Proxlike!*')

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
	if (wantsum == True):
		BigStopFile=open(WhichDir+'/bigstopfile.txt','w')
		BigStopFile.write(str(bigstop)+'\n')
		BigStopFile.close()

### Write summary to file
	done=False
#	if "Prx" not in WhichDir:
	if (wantsum == True):
		done=AC.WriteSummary(WhichDir,summary, summaryheader, stop, bigstop)

############################################################################
### Write summary to file
def WriteSummary(WhichDir, summary, summaryheader, stop, bigstop):

#	print('	WriteSummary '+WhichDir)

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
def MakePlots(version, WhichDir, t, eps, k, u, r, m, suffix=''):
	'''Make plots of parameters from this run'''

	print('	-MakePlots     '+WhichDir)

	mu = G*sum(m)

# Assign colors
	if version=='binary':
		c=('b','y')
	elif version=='triple':
		c=('y','r')
	else:
		print('MakePlots version name invalid! Pick "binary" or "triple".')

# Plot energies vs. separation
	plt.plot(r/AU, eps, 'ks')
	plt.plot(r/AU,   u, 'bo',
			 r/AU,   k, 'r^',
			)
	plt.legend(('Epsilon total','Potential','Kinetic'),
	           'lower right')
	plt.xlabel('Separation (AU)')
	plt.ylabel('Specific energy (J/kg)')
	plt.title('Energies')
	plt.savefig(WhichDir+'/EpsvsR'+suffix+'.png')
	plt.clf()
	
# Plot energies over time
	plt.plot(t, eps, 'k-')
	plt.plot(t,   u, 'b-',
			 t,   k, 'r-',
			)
	plt.legend(('Epsilon total','Potential','Kinetic'),
	           'lower right')
	plt.xlabel('time (years)')
	plt.ylabel('Energy (J)')
	plt.title('Energies')
	plt.xscale('log')
	plt.savefig(WhichDir+'/EpsvsT'+suffix+'.png')
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

 	InfoFile=open(WhichDir+'/Out/info.out','r')
	InfoLen=AC.FileLength(WhichDir+'/Out/info.out')
	AllInfo=InfoFile.readlines()
	InfoFile.close()

	start=where(AllInfo,"   Beginning the main integration.\n")
	ends=where(AllInfo,"   Integration complete.\n")
	nloops=len(ends)

	InfoBody = np.array(AllInfo[ (start[0]): ])

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
	FindTimeInd = []
	CompleteInd = []
	for i,row in enumerate(InfoBody):
		if (row != '\n'):
			if not ( (row.split()[0] == 'Fractional' ) | 
					 (row.split()[0] == 'Integration') | 
					 (row.split()[0] == 'Continuing')  | 
					 (row.split()[0] == 'WARNING:')    |
					 (row.split()[0] == 'Modify')      |
					 (row.split()[0] == 'Beginning' )   ):
				want.append(row)
			if (row.split()[0] == 'Continuing'):
				FindTimeInd.append(i)
			if (row.split()[0] == 'Integration' ):
				CompleteInd.append(i)

	# Find the start time of the last iteration
	FindTime = [float(InfoBody[i].split()[6]) for i in FindTimeInd]
	if len(CompleteInd)==0:
		if len(FindTime)>0:
			PrevTime = np.log10(max(FindTime))
		else:
			PrevTime = 0.
		complete = False
	elif(len(CompleteInd) == 1):
		PrevTime = 0.
		complete = True
#		print('First iteration, no time written'.format(PrevTime))
	elif (max(CompleteInd) > max(FindTimeInd)):
		PrevTime = max(np.log10(FindTime))
		complete = True
#		print('Last iteration complete, final time = {0}'.format(PrevTime))
	else:
		PrevTime = max(np.log10(FindTime))
		complete = False
#		print('Last iteration incomplete, final time = {0}'.format(PrevTime))

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

	

	return name,dest,time,10.**PrevTime,complete

###############################################################################
def SumAll(WhichDirs,cent,suffix=''):
	'''Read data from all directories and put in one file'''
	
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
#	sumpath=WhichDir+'/summary.out'
#
#	SumFile=open(sumpath, 'a')
#	if os.path.getsize(sumpath)==0:
#		SumFile.write('    aB     eB     iB      aC     eC     iC'+\
#		'         rBf          EBf         rCf          ECf'+\
#		'            tB    destB            tC    destC\n')
#	SumFile.close()
###############################################################################
def ReadDisk(WhichDir):
	'''Read in full time-dependent data on all objects from an AC Disk simulation'''
		
	print('        --ReadDisk     '+WhichDir)

	### Get list initial objects
	objsB = Merc.ReadInObjList(WhichDir+'/In/','big.in')
	center = objsB[0]   # this is the object the small bodies are orbiting
	objsS = Merc.ReadInObjList(WhichDir+'/In/','small.in',centerobj=center)
	namelist = np.array([o.name for o in objsB]+[o.name for o in objsS])
	alist    = [o.a()  for o in objsS]
	elist    = [o.e()  for o in objsS]
	mlist    = [o.mass for o in objsS]
	# Read in AEI data on disk objects
	nobjs  = len(namelist)  # big+small
	nparam = 6              # x y z vx vy vz
	ntime,time  = AC.GetNTime(WhichDir,namelist=namelist)     # n timesteps in longest aei file
	xv = np.empty((nobjs,ntime,nparam))
	xv.fill(np.nan)
	for ind, obj in enumerate(namelist):
		thisxv      = AC.ReadAei(WhichDir, obj, index1=None)
		xv[ind,0:thisxv.shape[0],:] = thisxv

	return(xv, namelist, time)

###############################################################################
def CalcDisk(WhichDir,center='default'):

	# Read in xv with dims (nobjs, ntime, nparam)
	xv, namelist, time = AC.ReadDisk(WhichDir)
	starmask = np.array([n     in ['AlCenA','AlCenB','PrxCen'] for n in namelist])
	diskmask = np.array([n not in ['AlCenA','AlCenB','PrxCen'] for n in namelist])
	diskxv = xv[diskmask,:,:]

	nobjs, ntimes, nparams = diskxv.shape
	# Calculate distance from central body for each obj at each time
	if center=='default':
		if   (center=='AlCenA') | (('A' in WhichDir) & ('B' not in WhichDir)):
			CentPos = np.zeros(xv[0,:,:].shape)
		elif (center=='AlCenB') | (('A' not in WhichDir) & ('B' in WhichDir)):
			CentPos = xv[0,:,:]
		else:
			assert 1==0, 'Error: uncertain center object; please specify "AlCenA" or "AlCenB"'

	### Binary stars' separation
#	rbin = np.sqrt( sum( [ xv[0,:,i]**2 for i in range(3) ] ) )
	### Objects' distance from central star
	r = np.sqrt( sum( [ (diskxv[:,:,i]-CentPos[:,i])**2 for i in range(3) ] ) )

	# Arrays for which objects are surviving/stable at each timestep
	surv = np.zeros((nobjs,ntimes),dtype='int')
	stab = np.zeros((nobjs,ntimes),dtype='int')
	for i in range(nobjs):
		for j in range(ntimes):
			if ( abs(r[i,j]-r[i,0])/r[i,0] < 0.15):
				surv[i,j] = 1
				stab[i,j] = 1
			elif ( not np.isnan(r[i,j]) ):
				surv[i,j] = 1
				stab[i,j] = 0
			else:
				surv[i,j] = 0
				stab[i,j] = 0

	### Calculate truncation radius over time
	# get index of last particle in disk at each timestep
	edge = np.array([AC.CalcDiskEdge(stab[:,i]) for i in range(ntimes)])

	# get truncation radius (original orbit of outermost particle in stable region
	# add last row of zeros, for cases with no disk left (edge=-1)
	rbound = np.concatenate(( r, np.zeros((1,ntimes)) ), axis=0 )
	rtr = np.array([rbound[edge[t],0] for t in range(ntimes)])

	return(WhichDir,time,r,rtr)

###############################################################################
### Commenting this function in favor of the next one.
### Doing these calculations with the whole array is more computationally
### efficient, but less intuitive, and I don't want to have multiple versions
### that could accidentally diverge with future editing
#def CalcDiskEdge2D(x, makeplot=False):
#	'''Takes an array of 0's and 1's and tries to fit a model in which the inner
#	part of the disk is mostly filled and the outer part is mostly empty. Returns
#	the index of the array which makes the best dividing point between them.'''

#	x = np.array(x)
#	nobjs  = x.shape[0]
#	ntimes = x.shape[1]
#	# Calculate avg stability on either side of each obj at each timestep
#	abov = np.zeros((nobjs,ntimes),dtype='float')
#	belw = np.zeros((nobjs,ntimes),dtype='float')
#	for i in range(nobjs):
#		belw[i,:] = np.mean(x[0:(i+1),:],axis=0)
#		if (i < nobjs-1):
#			abov[i,:] = np.mean(x[(i+1):,:],axis=0)
#		else:
#			abov[i,:] = 0.
#	# how closely does each point fit a solid disk with nothing beyond?
#	error = (1.-belw) + (abov)
#	# Fit edge of surviving disk
#	edge = np.zeros((ntimes),dtype='int')
#	for j in range(ntimes):
#		if all( error[:,j]==1 ): 
#			if all(belw[:,j]==0) and all(abov[:,j]==0):
#				if x[0,j] == 0:
#					# no disk left
#					edge[j] = -1
#				elif x[0,j] == 1:
#					# only innermost particle
#					edge[j] = 0
#			elif all(belw[:,j]==1) and all(abov[:,j]==1):
#				# full disk
#				edge[j] = nobjs-1 
#			else:
#				assert 0==1, 'Edge fitting algorithm is broken!'
#		else:
#			mininds = Merc.which( error[:,j],min(error[:,j]) )
#			if len(mininds)==1:
#				edge[j] = mininds[0]
#			elif len(mininds)>1:
#				# if multiple minimums, take 5-point means around each and use lower
#				nearmin = np.array([np.mean( error[(ind-2):(ind+3),j] ) for ind in mininds])
#				reallyMinInds = np.array(Merc.which( nearmin, min(nearmin) ))
#				if len(reallyMinInds) == 1:
#					edge[j] = mininds[reallyMinInds[0]]
#				elif len(reallyMinInds)>1:
#					print('Finding Rtr: multiple minimums in row {0}, inds {1}'.format(j,reallyMinInds))
#					edge[j] = np.mean( mininds )
#	return(edge)
#	
###############################################################################
def CalcDiskEdge(x, makeplot=False):
	'''Takes an array of 0's and 1's and tries to fit a model in which the inner
	part of the disk is mostly filled and the outer part is mostly empty. Returns
	the index of the array which makes the best dividing point between them.'''

	x = np.array(x)
	assert len(x.shape)==1,'Array with more than one dimension supplied! x.shape = {}'.format(x.shape)
	nObjs = len(x)    # count initial objects
	nTot  = np.sum(x) # count surviving objects
	### Calculate truncation radius over time
	# Calculate avg stability on either side of each obj at each timestep
	abov = np.zeros(nObjs,dtype='float')
	belw = np.zeros(nObjs,dtype='float')
	for i in range(nObjs):
		belw[i] = np.mean(x[0:(i+1)])
		if (i < nObjs-1):
			abov[i] = np.mean(x[(i+1):])
		else:
			abov[i] = 0.
	### Find fraction of filled slots above this point compared to below
	nAbv = np.array([np.sum(x[(i+1): ]) for i in range(nObjs)], dtype='float')
	nBlw = np.array([np.sum(x[0:(i+1)]) for i in range(nObjs)], dtype='float')
	ratio = nAbv/nBlw
	# how closely does each point fit a solid disk with nothing beyond?
	error = (1.-belw) + (abov) + nAbv/nTot
	# Fit edge of surviving disk
	if nTot > 1:
		# find lowest point(s) in error f'n
		mininds = Merc.which( error ,min(error) )
		if len(mininds)==1:
			edge = mininds[0]
		elif len(mininds)>1:
			# if multiple minimums, take 5-point means around each and use lower
			nearmin = np.array([np.mean( error[(ind-2):(ind+3)] ) for ind in mininds])
			reallyMinInds = np.array(Merc.which( nearmin, min(nearmin) ))
			if len(reallyMinInds) == 1:
				edge = mininds[reallyMinInds[0]]
			elif len(reallyMinInds)>1:
				print('AC.CalcDiskEdge: multiple minimums in row {0}, inds {1}'.format(j,reallyMinInds))
				edge = np.mean( mininds )
	elif nTot == 1:
		# only one surviving object
		if x[0] == 1:
			# if in lowest orbit, it counts as a disk
			edge = 0
		else:
			# if somewhere higher, ignore it (this is to override the too-high 
			# weight from ratio when there's only one object in the denominator
			edge = -1
	else:
		# no surviving objects
		edge = -1
#		if all(belw==0) and all(abov==0):
#			if x[0] == 0:
#				# no disk left
#				edge = -1
#			elif x[0] == 1:
#				# only innermost particle
#				edge = 0
#		elif all(belw==1) and all(abov==1):
#			# full disk
#			edge = nObjs-1 
#		else:
#			assert 0==1, 'Edge fitting algorithm is broken!'


	### Plot error f'ns and selected point, if requested
	if (makeplot==True):
		print(error)
		plt.clf()
		plt.plot(    error,'k',  label='error')
		plt.plot(     abov,'g',  label='abov')
		plt.plot(  1.-belw,'r',  label='belw')
		plt.plot(nAbv/nTot,'b',  label='ratio')
		plt.plot(edge, error[edge], 'ko', ms=10, label='edge')
		plt.legend()

	return(edge)

###############################################################################
def PlotDisk((WhichDir,time,r,rtr), samprate=1000, edgeline=True):
	'''Plot Disk data. First input is (CalcDisk()) output.'''

	nobjs, ntimes = r.shape
	
	### Assign colors for plot based on each object's stability over time
	# ndisk numbers spaced linearly from 0 to 1
	colspacing = linspace(1./3., 1., nobjs)
	# use above spacing to generate ndnobjsisk colors from heat colormap
	disk_cols = np.array([ cm.gist_heat(x) for x in colspacing ])
	# rgba of needed colors
	star_col  = colors.hex2color(colors.cnames['blue'])+(1.,)
	unst_col  = colors.hex2color(colors.cnames['lightgray'])+(1.,)
	ejec_col  = colors.hex2color(colors.cnames['gray'])+(1.,)

	# fill an array with the color for each object/time
	cols = np.empty((nobjs,ntimes,4.))
	for i in range(nobjs):
		for j in range(ntimes):
			if ( abs(r[i,j]-r[i,0])/r[i,0] < 0.15):
				cols[i,j,:] = disk_cols[i]
			elif ( not np.isnan(r[i,j]) ):
				cols[i,j,:] = unst_col
			else:
				cols[i,j,:] = ejec_col

	# Resample with log10 spacing to plot:
	# set up time spacing
	logt_goal = log10(time[1]) + \
			(log10(time[-1])-log10(time[1]))* \
			np.array(range(samprate)).astype(float)/(samprate-1.)
	t_goal = 10.**logt_goal
	# find indices in old time array that closest match resampled times
	tind_resamp = np.zeros((samprate),dtype='int')
	for j in range(samprate):
		tdist   = np.abs(t_goal[j]-time)
		tind_resamp[j] = np.array(range(len(time)),dtype='int')[tdist == min(tdist)][0]
	t_resamp = time[tind_resamp]
	logt_resamp = np.log10(t_resamp)
	# find values of other parameters at the resampled times
	rtr_resamp = rtr[tind_resamp]
	# Resampled image colors and rtr
	img = np.empty( (nobjs, samprate, 4) )
	for i in range(nobjs):
		for j,t_ind in enumerate(tind_resamp):
			# find closest value of time to resampt
			img[i,j,:] = cols[i, t_ind, :]

	# Set up figure
	f, ax1 = plt.subplots(2,figsize=(7,9))
	# grid of colors
	imglims = [ min(logt_resamp), max(logt_resamp), min(r[:,0]), max(r[:,0]) ]
	ax1[0].imshow(img[:,:], interpolation='nearest', aspect='auto',
		origin='lower', extent = imglims)
	ax1[0].set_xlabel('log10(t) (yrs)')
	ax1[0].set_ylabel('r/p_bin')	# pbin changes over time -- is this final?
	ax1[0].set_xlim(imglims[:2])
	ax1[0].set_ylim(imglims[2:])

	# trace truncation radius
	if edgeline==True:
		ax1[0].plot(logt_resamp,rtr_resamp)

	# colored lines
	for i,ri in enumerate(r):
		ax1[1].plot(time, ri, color=cols[i,-1,:])

	ax1[1].set_xscale('log')
#	ax1[1].set_yscale('log')

	ax1[1].set_ylim((0.1,11))

	ax1[1].set_xlabel('t')
	ax1[1].set_ylabel('a')

#	plt.show()
	f.savefig(WhichDir+'/Disk.png')
	f.clf()
	plt.close(f)

###############################################################################
def GetNTime(WhichDir,AeiDir='/Out/AeiOutFiles/',namelist='search'):
	'''Get number of timesteps and the time array from the longest .aei file'''

	AeiDir = WhichDir+AeiDir
	### If namelist not specified, get list of files in AeiOutDir
	if np.all(namelist=='search'):
		for root, dirs, files in os.walk(AeiDir):
			filelist=files
		filelist=np.array(filelist)
		validname = np.array([('.aei' in i) for i in filelist])
		filelist = filelist[validname]
	else:
		filelist = np.array([n+'.aei' for n in namelist])

### Get list of file sizes
	size=[]
	for f in filelist:
		size.append(os.path.getsize(AeiDir+f))
	size=np.array(size)
	### Get table dimensions of first on the list of max-length .aei files
	longfiles = filelist[size == np.max(size)]
	xv = AC.ReadAei(WhichDir, longfiles[0].strip('.aei'), None)
	ntime=xv.shape[0]
	# Get time array from this file
	time = AC.GetT(WhichDir,longfiles[0].strip('.aei'), 4,None)

	return(ntime,time)

###############################################################################





