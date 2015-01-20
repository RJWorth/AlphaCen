###############################################################################
def ComparePipedParams(WhichDir,cases=['O','A','B']):
	'''Compare orbital parameters'''
	
### Modules
	import AlphaCenModule as AC
	import numpy as np
	import Compare

### Currently, must have original version included
	assert 'O' in cases

	print(WhichDir)
### Read in the three relevant run.pipe files
	nOlines = AC.FileLength(WhichDir+'/Original/run.pipe')
	origfile=open(          WhichDir+'/Original/run.pipe')
	O=origfile.readlines()
	origfile.close()

	if 'A' in cases:
		nAlines = AC.FileLength(WhichDir+'/DiskA-3/run.pipe')
		Afile=open(          WhichDir+'/DiskA-3/run.pipe')
		A=Afile.readlines()
		Afile.close()
		AaeiB,AaeiC,AT = Compare.FindParams(A,nAlines)
	else:
		AT = 0.

	if 'B' in cases:
		nBlines = AC.FileLength(WhichDir+'/DiskB-3/run.pipe')
		Bfile=open(          WhichDir+'/DiskB-3/run.pipe')
		B=Bfile.readlines()
		Bfile.close()
		BaeiB,BaeiC,BT = Compare.FindParams(B,nBlines)
	else:
		BT = 0.

### Check lengths against each other
	if ('A' in cases) & ('B' in cases):
		if (nAlines != nBlines):
			print('Warning: lengths of run.pipe for A-3 and B-3 do not match!')
### Find param lines from end of A
		if (AT != BT):
			print('Warning: lengths of run.pipe for A-3 and B-3 do not match! {0} {1}'.format(AT,BT))
	if ( max(float(AT),float(BT)) == float(AT) ):
		T = AT
	else:
		if ('A' in cases) & ('B' not in cases):
			T = AT
		else:
			T = BT
	if ('071714' in WhichDir):
		print('071714: using last timestep from Original')
		OaeiB,OaeiC,OT = Compare.FindParams(O,nOlines)	
	else:
		OaeiB,OaeiC,OT = Compare.FindParams(O,nOlines,T)

	allT = np.array([])
	allT = np.append(allT,OT)
	Binary = [OaeiB]
	Triple = [OaeiC]
	if 'A' in cases:
		allT=np.append(allT,AT)
                Binary=np.vstack((Binary,AaeiB))
                Triple=np.vstack((Triple,AaeiC))
	if 'B' in cases:
		allT=np.append(allT,BT)
                Binary=np.vstack((Binary,BaeiB))
                Triple=np.vstack((Triple,BaeiC))

	display=np.c_[cases,Binary,Triple,allT]
	print(['','   aB','    eB','   iB','   aC','   eC','   iC','t'])
	print(display)

	names=['a','e','i']
	for i in range(len(cases)):
		if (allT[i]!=allT[0]):
			print('Mismatch in time: {0}'.format(allT))
	matching=6
	for j in range(3):
		for i in range(1,len(cases)):
			a,b = float(Binary[i,j].strip(',')), float(Binary[0,j].strip(','))
			eps = abs(a-b)/abs(min(a,b))
#			print(j,i,a,b,eps)
			if ( (eps > 0.01) | np.isnan(eps) ):
				print('***{0:.1f}% mismatch in {1}: {2}***'.format(eps*100.,names[j],Binary[:,j]))
				matching=matching-1
			ab,c= float(Triple[i,j].strip(',')), float(Triple[0,j].strip(','))
			eps2= abs(ab-c)/abs(min(ab,c))
#			print(j,i,ab,c,eps2)
			if ((eps2 > 0.01) | np.isnan(eps2)):
				print('***{0:.1f}% mismatch in {1}: {2}***'.format(eps2*100.,names[j],Triple[:,j]))
				matching=matching-1

### Is the new sim also proxlike?
	apo = float(BaeiC[0].strip(',')) * ( 1 + float(BaeiC[1].strip(',')) )
	if (apo>10000):
		proxlike=1
	else:
		proxlike=0 
	print('Proxlike = ',proxlike)
	
	Matchfile=open(WhichDir+'/match.txt','w')
	Matchfile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}".format(
		matching,proxlike,OT,
		OaeiB[0].strip(','),OaeiB[1].strip(','),OaeiB[2].strip(','),
		OaeiC[0].strip(','),OaeiC[1].strip(','),OaeiC[2].strip(',') ) )
#	Matchfile.write(str(matching))	
	Matchfile.close()
	

###############################################################################
def FindParams(A,nlines,time='default'):
	'''Get B and C orbital parameters and end time from run.pipe'''

	import numpy as np

	stop1,stop2,stop3 = False,False,False
	for i in range(nlines):
		j=nlines-i-1
		if ( len(A[j].split()) > 0 ):
			if (A[j].split()[0] == 'stop'):
				if (time=='default'):
					indT=j
					stop3=True
				elif (A[j].split()[5] == time):
					indT=j
					stop3=True
			if (stop3 == True):
				if (  A[j].split()[0] == 'aC'):
					indC=j
					stop1=True
				elif (A[j].split()[0] == 'aB'):
					indB=j
					stop2=True
			if (stop1 & stop2 & stop3):
				break
			if (i==nlines-1):
				print('Incomplete matching, stops = {0} {1} {2}'.format(
						stop1,stop2,stop3))
	if (stop1 & stop2 & stop3):
		aeiB = np.array(A[indB].split())[ [2,5,8] ]
		aeiC = np.array(A[indC].split())[ [2,5,8] ]
		T = A[indT].split()[ 5 ]
	else:
		aeiB,aeiC,T = ['nan', 'nan','nan'], [ 'nan','nan','nan' ], 0.

	return (aeiB,aeiC,T)
