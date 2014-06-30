###############################################################################
### Function to count the number of lines in a file
def FileLength(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

###############################################################################
# Read .aei file to get xyz and uvw (in AU and m/s respectively)
def ReadAei(whichdir, filename, index1=-1, index2=0):

	print('ReadAei '+whichdir+'/Out/AeiOutFiles/'+filename+', '+\
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
# Read .aei file to get xyz and uvw (in AU and m/s respectively)
def GetT(whichdir, filename, index1=-1, index2=0):

	print('ReadAei '+whichdir+'/Out/AeiOutFiles/'+filename+', '+\
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
# Get distance between two points
def Distance(xv1, xv2):

	r = (  (xv1[:,0]-xv2[:,0])**2
		 + (xv1[:,1]-xv2[:,1])**2
		 + (xv1[:,2]-xv2[:,2])**2	)**0.5

	return(r)

###############################################################################
# Take xv vector(s), return r
def XVtoR(xv):

	r = ( xv[:,0]**2 + xv[:,1]**2 + xv[:,2]**2 )**0.5

	return(r)

###############################################################################
# Take xv vector(s), return v
def XVtoV(xv):

	v = ( xv[:,3]**2 + xv[:,4]**2 + xv[:,5]**2 )**0.5

	return(v)

###############################################################################
### Estimate a from v(r) function
def RVtoA(r, v, mu):
# where r is the distance between two stars, v is the relative velocity

	a = 1/( 2/r -  v**2/mu )

	return(a)

###############################################################################
### Estimate v from v(r) function
def RAtoV(r, a, mu):
# where r is the distance between two stars, v is the relative velocity

	v = (mu*( 2/r -  1/a) )**0.5

	return(v)

###############################################################################
# Convert xv in AU, AU/day to mks
def AUtoMKS(xv_AU):

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
### Take position and velocity data, give center of mass/momentum versions
def FindCM(m, xv):

	print('FindCM: m = '+str(m))

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
### Transform given coordinates to CM frame
def wrtCM(xv, xvCM):

	print('wrtCM')

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
### Calculate kinetic energy
def Kinetic(m, v):
# where r is the distance between two stars, v is the relative velocity
	import numpy as np

	K = [(1./2.)*m[i]*v[i,:]**2. for i in range(len(v[:,0]))]

	return(np.array(K))

###############################################################################
### Calculate potential energy of one object due to another
def Potential(m1, m2, r):
# where r is the distance between two stars, v is the relative velocity
	from mks_constants import G

	U = -G*m1*m2/r
	U = U/2.

	return(U)

###############################################################################

