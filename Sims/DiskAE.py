###############################################################################
def WriteAEI(WhichDir, n):

#WhichDir = 'Proxlike/Prx01/DiskB-2'
#n = 10

	import AlphaCenModule as AC
	import numpy as np
	from mks_constants import G, mSun, AU, day

### Get object list for stars + n M-particles
	if 'B-3' in WhichDir:
		SmlInd=3
		objs = ['AlCenA','AlCenB','PrxCen']+['M'+str(i) for i in range(1,n+1)]
		if (any([s in WhichDir for s in ['071714','072314','081414']] ) ):
			m = np.array([1.105,0.934,0.123])*mSun
		else:
			m = np.array([0.123,0.123,0.123])*mSun		
	elif 'B-2' in WhichDir:
		SmlInd=2
		objs = ['AlCenA','AlCenB']+['M'+str(i) for i in range(1,n+1)]
		if (any([s in WhichDir for s in ['071714','072314','081414']] ) ):
			m = np.array([1.105,0.934])*mSun
		else:
			m = np.array([0.123,0.123])*mSun		

### Get corresponding time array and maximum # of timesteps from AlCenB
	t = AC.GetT(WhichDir, objs[1], 4, 0)
	print(t)
	np.savetxt(WhichDir+'/Out/AeiOutFiles/t.out',np.transpose(t))
	maxT = len(t)

### Read in AlCenB
	xvB_AU = AC.ReadAei(WhichDir, objs[1], None,None)
	xvB = AC.AUtoMKS(xvB_AU)

### Make 3D arrays to fill, based on these measurements
	xv  = np.zeros( (len(objs), maxT, 6) )
	aei = np.zeros( (len(objs), maxT, 3) )
	a = np.zeros( (len(objs), maxT) )
	e = np.zeros( (len(objs), maxT) )
	i = np.zeros( (len(objs), maxT) )
	x = np.zeros( (len(objs), maxT) )
	y = np.zeros( (len(objs), maxT) )
	z = np.zeros( (len(objs), maxT) )

### Fill AlCenB slice
	xv[1,:,:] = xvB

### Fill remainder of array
	for j in range(2,len(objs)):
		thisxv_AU = AC.ReadAei(WhichDir, objs[j], None,None)
		thisT     = thisxv_AU.shape[0]
		thisxv    = AC.AUtoMKS(thisxv_AU)
		xv[j,:thisT,:] = thisxv
	
### Make arrays for x, y, z components
	x=xv[:,:,0]
	y=xv[:,:,1]
	z=xv[:,:,2]

### Get AB binary orbit
	aB, eB, iB, epsB, xvCM_AB, kB, uB, rB, vB = AC.Binary(
				[ m[0], m[1] ], 
				np.array([ xv[0,:,:], xv[1,:,:] ]), 
				maxT)
	aei[1,:,0],aei[1,:,1],aei[1,:,2] = aB/AU, eB, iB
	a[1,:]    ,    e[1,:],    i[1,:] = aB/AU, eB, iB

### Get Prx orbit relative to CM_AB (if applicable)
	if 'B-3' in WhichDir:
		aP, eP, iP, eps, xvCM, k, u, r, v = AC.Binary(
				[ m[0]+m[1], m[2] ], 
				np.array([ xvCM_AB, xv[2,:,:] ]), 
				maxT)
		aei[2,:,0],aei[2,:,1],aei[2,:,2] = aP/AU, eP, iP
		a[2,:]    ,    e[2,:],    i[2,:] = aP/AU, eP, iP

### Get small obj orbits relative to B
	for j in range(SmlInd,len(objs)):
		aM, eM, iM, eps, xvCM, k, u, r, v = AC.Binary(
				[ m[1], 0. ], 
				np.array([ xv[1,:,:], xv[j,:,:] ]), 
				maxT)
		aei[j,:,0],aei[j,:,1],aei[j,:,2] = aM/AU, eM, iM
		a[j,:]    ,    e[j,:],    i[j,:] = aM/AU, eM, iM

	np.savetxt(WhichDir+'/Out/AeiOutFiles/x.out',np.transpose(x),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')
	np.savetxt(WhichDir+'/Out/AeiOutFiles/y.out',np.transpose(y),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')
	np.savetxt(WhichDir+'/Out/AeiOutFiles/z.out',np.transpose(z),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')

	np.savetxt(WhichDir+'/Out/AeiOutFiles/a.out',np.transpose(a),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')
	np.savetxt(WhichDir+'/Out/AeiOutFiles/e.out',np.transpose(e),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')
	np.savetxt(WhichDir+'/Out/AeiOutFiles/i.out',np.transpose(i),fmt='%.2e',
		header=''.join([o.rjust(9) for o in objs]),comments='')
	
	print(xv[:,-1,:])

###############################################################################
