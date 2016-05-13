###############################################################################
def NCA(a=23.5,e=0.52,mu=.46,Mstar=.93):
	'''Calculate N_CA based on Jang-Condell 2015 (Default for AlCenB)'''

	import numpy as np

### Make array of coefficients
	cCA=np.zeros((2,2,3,3))
	cCA[0,0,:,0] = [-6.7585, -2.0956,  1.8452]
	cCA[1,0,:,0] = [ 0.4308,  1.6243, -0.5589]
	cCA[0,1,:,0] = [ 2.9369, -9.5050, 16.0528]
	cCA[1,1,:,0] = [-0.4411, -1.4578,  0.3879]
	
	cCA[0,0,:,1] = [ 1.5486, -9.3964,  7.9117]
	cCA[1,0,:,1] = [-0.1609,  0.2216, -0.3244]
	cCA[0,1,:,1] = [-1.9967,  3.5221, -5.4470]
	cCA[1,1,:,1] = [ 0.1696, -0.3374,  0.4971]

	cCA[0,0,:,2] = [-0.2951,  2.7651, -2.2214]
	cCA[1,0,:,2] = [ 0.0364, -0.0993,  0.0920]
	cCA[0,1,:,2] = [ 0.4372, -1.6003,  2.1346]
	cCA[1,1,:,2] = [-0.0363,  0.1275, -0.1406]

### Calculate N_CA
	NCA = 0.
	for i in range(2):
		for j in range(2):
			for k in range(3):
				for l in range(3):
					term = cCA[i,j,k,l] * a**i * e**j * mu**k * Mstar**l
					NCA = NCA + term

	return(NCA)

###############################################################################
def NDI(a=23.5,e=0.52,mu=.46,Mstar=.93):
	'''Calculate N_CA based on Jang-Condell 2015 (Default for AlCenB)'''

	import numpy as np

### Make array of coefficients
	cDI=np.zeros((2,2,2,2))
	cDI[0,0,:,0] = [ 0.3921, -0.1135]
	cDI[1,0,:,0] = [ 0.0521,  0.1106]
	cDI[0,1,:,0] = [-0.8329,  1.1926]
	cDI[1,1,:,0] = [-0.0431, -0.0975]
	
	cDI[0,0,:,1] = [-0.2704,  0.2944]
	cDI[1,0,:,1] = [-0.0066, -0.0296]
	cDI[0,1,:,1] = [ 0.3131, -0.7153]
	cDI[1,1,:,1] = [ 0.0012,  0.0327]

### Calculate N_CA
	NDI = 0.
	for i in range(2):
		for j in range(2):
			for k in range(2):
				for l in range(2):
					term = cDI[i,j,k,l] * a**i * e**j * mu**k * Mstar**l
					NDI = NDI + term

	return(NDI)

###############################################################################
def mDust(ro,ri=0.35,alpha=1.5,sigma0cgs=1700.):
	'''Calculate the mass of a dust disk with slope alpha, density sigma0.cgs
	in g/cm^2 at 1 AU, inner edge ri, and outer/truncation radius ro.'''
### Make sure the inputs are written as floats, not integers

### Two alpha value cases allowed
	assert ((alpha==1.5) | (alpha==1))
	assert (ro>= ri)

	import numpy as np
	from numpy import pi
	from cgs_constants import mEarth, AU

### Convert sigma to mSun/AU^2
	sigma0 = sigma0cgs*AU**2./mEarth
	if alpha==1:
		sigma0=0.303*sigma0			# based on equal-mass 36-AU disks -- 5/4/15

### Calculate mTot in annulus with specified density profile
	if alpha==1.5:
		mTot = 4.*pi*sigma0* (1)**1.5 *(ro**0.5 - ri**0.5)
	elif alpha==1:
		mTot = 2.*pi*sigma0* (1)      *(ro      - ri)

	return mTot

###############################################################################
def mDustTot(ro,ri=0.35, alpha=1.5,sigcoef=1.,rice=2.7):
	'''Calculate the mass of a dust disk using the Hayashi slope formula.'''

	import numpy as np
	import JangCondell as J
	from numpy import pi
	from cgs_constants import mEarth, mSun, AU

### Densities, breakpoints from Hayashi '81 paper
	rho = sigcoef*np.array([7.1, 30., 1700.])
	rin  = 0.35	# inner edge
	rice = rice	# ice line
	rout = 36.	# outer edge

### Inside ice line
	if (ri < rice):
		ri1 = max( min(ri, rice), rin)
		ro1 = min( max(ro,  rin), rice)
		mDust1 = J.mDust( ro1,ri1,alpha=alpha, sigma0cgs=rho[0] ) 
	else:
		mDust1 = 0.
		print('no inner disk')

### Outside ice line
	if (ro > rice):
		ri2 = max( ri, rice)
		ro2 = min( ro, rout)
		mDust2 = J.mDust( ro2,ri2,alpha=alpha, sigma0cgs=rho[1] ) 
	else:
		mDust2 = 0.
		print('no outer disk')

### Gas component
	rig = max( ri, rin)
	rog = min( ro, rout)
	mGas   = J.mDust( rog,rig,alpha=alpha, sigma0cgs=rho[2] ) 

### Total solid mass
#	mDust = mDust1+mDust2

	return [mDust1, mDust2, mGas]


###############################################################################
def DiskTable():
	'''Calculate dust mass in truncated disks for a suite of parameters'''

	import numpy as np
	import JangCondell as J
	from numpy import pi
	from cgs_constants import mEarth, mSun, AU


### Parameters
	rice    = np.array([ 2.5, 2.7, 3.])
	sigcoef = np.array([1/3.,  1., 3.])
	alpha   = np.array([  1., 1.5])
	rtr     = np.array([2.54, 2.77])

### Make 3x3x2x2 array, fill with mDust for above parameters
	ind = np.zeros( shape=( len(rtr), len(alpha), len(sigcoef), len(rice) ) )
	M   = np.zeros( shape=( len(rtr), len(alpha), len(sigcoef), len(rice) ) )
	counter=0
	for l in range(len(rtr)):
		for k in range(len(alpha)):
			for j in range(len(sigcoef)):
				for i in range(len(rice)):
					counter=counter+1
					ind[l,k,j,i] = counter
					mTot = J.mDustTot(
				ro=rtr[l], alpha=alpha[k],sigcoef=sigcoef[j],rice=rice[i])
					M[l,k,j,i] = mTot[0]+mTot[1]

### Collapse over two dimensions
	M = np.hstack(( M[:,0,:,:], M[:,1,:,:] )) 
	M = np.hstack(( M[0,:,:],   M[1,:,:]   )) 

	np.savetxt("../Paper/Inserts/DiskMassTable.tex", M, 
		fmt = '%2.2f', delimiter=' & ', newline=' \\\\\n')

	return(M)




