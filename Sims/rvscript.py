
### Modules
import rvtest as rv
from mks_constants import G, mSun, AU, day, m, mu
import numpy as np
from numpy import sin, pi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# mass parameters
m = np.array(m)
m2 = m[0:2]			# binary only
m3 = m[0:3]			# triple
mu2=sum(m[0:2])*G	# binary only
mu3=sum(m[0:3])*G	# triple
kfactor=G*day**2/AU**3	# kfactor*m = k^2 in AU^3/d^2

### Import controls
WhichDir='TestD/Out/AeiOutFiles/'
suffix=''
index1=-1000
index2=0
mode='binary'	# 'binary' or 'triple'

###
assert (mode=='triple') | (mode=='binary')

if (mode == 'binary'):
	m=[m[0], m[1], 0.]
	filenames=['AlCenB']
elif (mode == 'triple'):
	filenames=['AlCenB', 'PrxCen']

nobjs=len(filenames)+1
ntB=rv.FileLength(WhichDir+filenames[0]+suffix+'.aei')-4
if (mode == 'triple'):
	ntC=rv.FileLength(WhichDir+filenames[1]+suffix+'.aei')-4

### NOTATION:
# xv = the position-velocity vector (x, y, z, u, v, w),
# or an array of these vectors over time, where
# xv[:,0] = the x column and
# xv[i,:] = xv at timestep i

# 3D array: xvA = xv of all objects w.r.t. CMA (in AU units)
# xvA[i,j,k] = xv of object i, time j, column k

# Units are mks unless specified otherwise by tacking _AU on the end, 
# then it's in the AU units that MERCURY outputs (x=AU, v=AU/day).

# rAB = distance between A and B
# rCMAB = distance of each star from the CM of A and B
# vCMAB = speed of each star, w.r.t. the CM of A and B

# xvI_J = xv of star I with respect to the center of mass of star(s) J

# get output times
t = rv.GetT(WhichDir, filenames[0]+suffix, index1, index2)

############# Read in original x and v values ###########################
xvB_A_AU		= rv.ReadAei(WhichDir, filenames[0]+suffix, index1, index2)
xvA_A_AU		= np.zeros_like(xvB_A_AU)
if (mode == 'triple'):
	xvC_A_AU	= rv.ReadAei(WhichDir, filenames[1]+suffix, index1, index2)
### Combine three stars into a 3D array
	xvA_AU=np.array([xvA_A_AU, xvB_A_AU, xvC_A_AU])
else:
### Or two stars, if just a binary system
	xvA_AU=np.array([xvA_A_AU, xvB_A_AU])

### 1st dimension of array should be the number of stars
assert(np.shape(xvA_AU)[0]) == nobjs

### Get distances between stars
#rAB_AU = rv.Distance(xvA_AU[0,:,:], xvA_AU[1,:,:])
#if (mode=='triple'):
#	rBC_AU = rv.Distance(xvA_AU[1,:,:], xvA_AU[2,:,:])
#	rAC_AU = rv.Distance(xvA_AU[0,:,:], xvA_AU[2,:,:])

### Get r, v in AU units
#rCMA_AU = np.array([ rv.XVtoR(xvA_AU[i,:,:]) for i in range(nobjs) ])
#vCMA_AU = np.array([ rv.XVtoV(xvA_AU[i,:,:]) for i in range(nobjs) ])

### Estimate a from v(r) (wrt other star)
#aCMA_AU   = np.array([ rv.RVtoA(rAB_AU, vCMA_AU[1,:], mu2*day**2/AU**3),
#					   rv.RVtoA(rAB_AU, vCMA_AU[1,:], mu2*day**2/AU**3) ])
#if (mode=='triple'):
#	aCMA_AU = np.array([ aCMA_AU, 
#			  rv.RVtoA(rAC_AU, vCMA_AU[2,:], mu3*day**2/AU**3) ])

### Get expected velocity (wrt other star)
# should be same as vCMA_AU of B, for both, because it's relative to the other
#vexpCMA_AU	= rv.RAtoV(rAB_AU, aCMA_AU, mu2*day**2/AU**3)


##################### Convert to mks units ##############################
xvA = rv.AUtoMKS(xvA_AU)

### Get distances between stars (== rAij_AU*AU)
rAB = rv.Distance(xvA[0,:,:], xvA[1,:,:])
if (mode=='triple'):
	rBC = rv.Distance(xvA[1,:,:], xvA[2,:,:])
	rAC = rv.Distance(xvA[0,:,:], xvA[2,:,:])

### Get r, v in mks
#rCMA = np.array([ rv.XVtoR(xvA[i,:,:]) for i in range(nobjs) ])
#vCMA = np.array([ rv.XVtoV(xvA[i,:,:]) for i in range(nobjs) ])

### Estimate a from v(r)
#aCMA	= np.array([ rv.RVtoA(rAB, vCMA[1,:], mu2),
#				     rv.RVtoA(rAB, vCMA[1,:], mu2) ])
#if (mode=='triple'):
#	aCMA = np.array([ aCMA, rv.RVtoA(rAC, vCMAB[2,:], mu3) ])

### Get expected velocity
# should be same as vCMA of B, for both, because it's relative to the other
#vexpCMA	= rv.RAtoV(rAB, aCMA, mu2)

### Energy in A-centered coords?
#KCMA = rv.Kinetic(m,vCMA)
#UCMA = np.array([ rv.Potential(m[0],m[1], rAB),
#				  rv.Potential(m[0],m[1], rAB) ])
#ECMA = KCMA+UCMA
#EtotCMA = np.sum(ECMA, 0)

######################## Binary system: #################################
# CMAB = Center of momentum frame of A+B
### Find the coordinates of the center of momentum
xvCM_AB = rv.FindCM( m[0:2], xvA[0:2,:,:]) 
### Convert to center-of-momentum units
xvAB  = rv.wrtCM(xvA, xvCM_AB)
### Get r, v in CM units
rCMAB = np.array([ rv.XVtoR(xvAB[i,:,:]) for i in range(nobjs) ])
vCMAB = np.array([ rv.XVtoV(xvAB[i,:,:]) for i in range(nobjs) ])
######################## CM consistency check ###########################

#dxvA  =  xvA[0,:] -  xvA[1,:]
#dxvAB = xvAB[0,:] - xvAB[1,:]

#print(dxvA-dxvAB)	# == 0. for everything


########################## Energies #####################################
# Get kinetic energy = (1/2)mv^2
KCMAB = rv.Kinetic(m,vCMAB)
# Get potential energy = GMm/r
UCMAB = np.array([ rv.Potential(m[0],m[1], rAB),
				   rv.Potential(m[0],m[1], rAB) ])

# Get total energy = K+U
ECMAB = KCMAB+UCMAB

### This should be constant, but isn't...?
EtotCMAB = np.sum(ECMAB, 0)

########################## Energies in AU units #########################
#xvCM_AB_AU = rv.FindCM( m[0:2]*kfactor, xvA_AU[0:2,:,:]) 
#xvAB_AU  = rv.wrtCM(xvA_AU, xvCM_AB_AU)
#vCMAB_AU = np.array([ rv.XVtoV(xvAB_AU[i,:,:]) for i in range(nobjs) ])
# Get kinetic energy = (1/2)mv^2
#KCMAB_AU = rv.Kinetic(m*kfactor,vCMAB_AU)
# Get potential energy = GMm/r
#UCMAB_AU = np.array([ rv.Potential(m[0]*kfactor/G,m[1]*kfactor, rAB_AU),
#				   rv.Potential(m[0]*kfactor/G,m[1]*kfactor, rAB_AU) ])

# Get total energy = K+U
#ECMAB_AU = KCMAB_AU+UCMAB_AU/2

### This should be constant, but isn't...?
#EtotCMAB_AU = np.sum(ECMAB_AU, 0)

#E_array=np.array((UCMAB_AU[0,:],np.sum(KCMAB_AU,0),ECMAB_AU[0,]))
#compare=np.array((np.sum(KCMAB_AU,0),vCMAB_AU[0,:],vCMAB_AU[1,:] ))

###################### Triple system: ###################################
### Only relevant if B and C both survived
if (mode=='triple'):
	if (ntC==ntB):
### Find the coordinates of the center of momentum
		xvCM_ABC = rv.FindCM(m, xvA) 
### Convert to center-of-momentum units
		xvABC  = rv.wrtCM(xvA, xvCM_ABC)
### Get r, v in CM units
		rCMABC = np.array([ rv.XVtoR(xvABC[i,:,:]) for i in range(nobjs) ])
		vCMABC = np.array([ rv.XVtoV(xvABC[i,:,:]) for i in range(nobjs) ])
########################## Energies #####################################
### Calculate kinetic, potential, and total energies
		K = rv.Kinetic(m, vCMABC)
### Calculate a, e of orbit
		UAB = rv.Potential(m[0],m[1], rAB)
		UBC = rv.Potential(m[0],m[2], rAC)
		UAC = rv.Potential(m[1],m[2], rBC)
		U = np.array([ UAB+UAC, UAB+UBC, UBC+UAC ])
# Calculate total energy per object K+U
		E = K+U
# Calculate total energy in the system
		Etot = np.sum(E, 0)

################################# Plot ##################################
### Plot
# v(r) and expected v(r)
#plt.plot(rAB/AU,   vCMA[1,:], 'ko',
#		 rAB/AU,vexpCMA[1,:], 'b^',
#		 rCMAB[2,:]/AU,vCMAB[2,:], 'rs',
#		 )
#plt.legend(('Actual','Expected'))
#plt.xlabel('A-B separation (AU)')
#plt.ylabel('v (m/s)')
#plt.title('CM = A')
#plt.savefig(WhichDir+'VvsR_CMA'+suffix+'.png')
#plt.clf()

# v(r) in CM frame
plt.plot(rAB/AU,   vCMAB[0,:], 'b^',
		 rAB/AU,   vCMAB[1,:], 'y^',
#		 rCMAB[2,:]/AU,vCMAB[2,:], 'rs',
		 )
plt.legend(('A','B'))
plt.xlabel('A-B separation (AU)')
plt.ylabel('v (m/s)')
plt.title('CM = A+B')
plt.savefig(WhichDir+'VvsR_CMAB'+suffix+'.png')
plt.clf()

# 
if (mode=='triple'):
	plt.plot(rAB/AU,        Etot, 'ko',
			 rAB/AU,      E[0,:], 'b^',
			 rAB/AU,      E[1,:], 'y^',
			 rAB/AU,      E[2,:], 'r^',
			 )
	plt.savefig(WhichDir+'EvsR_Triple'+suffix+'.png')
	plt.clf()

# Plot energies of binary
rtest = np.array(range(19,37))*AU
normalization = np.array(range(-8,3))*1e37
Etest=np.empty( (len(normalization),len(rtest)) )
for i in range(len(normalization)):
	Etest[i,:] = normalization[i]/(rtest/np.max(rAB))

Ktest = np.min(KCMAB)/(rtest/np.max(rAB))
Utest = np.max(UCMAB)/(rtest/np.max(rAB))


plt.plot(rAB/AU,   EtotCMAB, 'ks',
		 rAB/AU, ECMAB[0,:], 'bs',
		 rAB/AU, ECMAB[1,:], 'ys',
		 rAB/AU, UCMAB[0,:], 'bo',
		 rAB/AU, KCMAB[0,:], 'b^',
		 rAB/AU, UCMAB[1,:], 'yo',
		 rAB/AU, KCMAB[1,:], 'y^',
		 rtest/AU, Etest[0,:],
		 rtest/AU, Etest[1,:],
		 rtest/AU, Etest[2,:],
		 rtest/AU, Etest[3,:],
		 rtest/AU, Etest[4,:],
		 rtest/AU, Etest[5,:],
		 rtest/AU, Etest[6,:],
		 rtest/AU, Etest[7,:],
		 rtest/AU, Etest[8,:],
		 rtest/AU, Etest[9,:],
		 rtest/AU, Etest[10,:],
		 )
plt.legend(('System total','A total','B total','A potential','A kinetic'),
           'lower right')
plt.xlabel('A-B separation (AU)')
plt.ylabel('Energy (J)')
plt.title('Energies, CM = A+B')
plt.savefig(WhichDir+'EvsR_CMAB'+suffix+'.png')
plt.clf()

# Plot energies over time
#t=1000*(np.array(range(1000))+1)
#moret=100*(np.array(range(10000))+1)
#sinfit = (-.21*sin((2*pi/6000.)*moret + 0.0) - 1.03)*1.e38
Amp		= .205
Freq	= 102.
Phase	= 0.3
Offset	= -1.05
sinfit = (Amp*sin((2*pi/Freq)*t + Phase) + Offset)*1.e38

plt.plot(t, sinfit, 
		 t, EtotCMAB, 'ks',
        )
#plt.xlim([0,21000])
#plt.legend(('System total'))
plt.xlabel('time (years)')
plt.ylabel('Energy (J)')
plt.title('Energies, CM = A+B')
plt.savefig(WhichDir+'EvsT_CMAB'+suffix+'.png')
plt.clf()

# Plot separation over time
plt.plot(t, rAB/AU, 'ks',
        )
#plt.xlim([0,21000])
#plt.legend(('System total'))
plt.xlabel('time (years)')
plt.ylabel('rAB (AU)')
plt.title('A-B separation')
plt.savefig(WhichDir+'RvsT_CMAB'+suffix+'.png')
plt.clf()


# Plot energies of binary in A-centered coords
#plt.plot(rAB/AU,   EtotCMA, 'ks',
#		 rAB/AU, ECMA[0,:], 'bs',
#		 rAB/AU, ECMA[1,:], 'ys',
#		 rAB/AU, UCMA[0,:], 'bo',
#		 rAB/AU, KCMA[0,:], 'b^',
#		 rAB/AU, UCMA[1,:], 'yo',
#		 rAB/AU, KCMA[1,:], 'y^',
#		 )
#plt.legend(('System total','A total','B total','A potential','A kinetic'))
#plt.xlabel('A-B separation (AU)')
#plt.ylabel('Energy (J)')
#plt.title('Energies, CM = A')
#plt.savefig(WhichDir+'EvsR_CMA'+suffix+'.png')
#plt.clf()

