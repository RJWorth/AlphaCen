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

		print('*** READING B-2 ***')
		nB2lines = AC.FileLength(WhichDir+'/DiskB-2/run.pipe')
                B2file=open(          WhichDir+'/DiskB-2/run.pipe')
                B2=B2file.readlines()
                B2file.close()
                B2aeiB,B2aeiC,B2T = Compare.FindParams(B2,nB2lines)
		print('*** DONE W/ B-2 ***')
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
			if ( (eps > 0.01) | np.isnan(eps) ):
				print('***{0:.1f}% mismatch in {1}: {2}***'.format(eps*100.,names[j],Binary[:,j]))
				matching=matching-1
			ab,c= float(Triple[i,j].strip(',')), float(Triple[0,j].strip(','))
			eps2= abs(ab-c)/abs(min(ab,c))
			if ((eps2 > 0.01) | np.isnan(eps2)):
				print('***{0:.1f}% mismatch in {1}: {2}***'.format(eps2*100.,names[j],Triple[:,j]))
				matching=matching-1
### Calculate change in a and e from Original to B-3 sim
	da = float(B2aeiB[0].strip(',')) - float(BaeiB[0].strip(','))
	de = float(B2aeiB[1].strip(',')) - float(BaeiB[1].strip(','))

### Is the new sim also proxlike?
	apo = float(BaeiC[0].strip(',')) * ( 1 + float(BaeiC[1].strip(',')) )
	if (apo>10000):
		proxlike=1
	else:
		proxlike=0 
	print('Proxlike = ',proxlike)
	
	Matchfile=open(WhichDir+'/match.txt','w')
	Matchfile.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}".format(
		matching,proxlike,OT,
		BaeiB[0].strip(','),BaeiB[1].strip(','),BaeiB[2].strip(','),
		BaeiC[0].strip(','),BaeiC[1].strip(','),BaeiC[2].strip(','),
		da, de ) )
	Matchfile.close()
	

###############################################################################
def FindParams(A,nlines,time='default'):
	'''Get B and C orbital parameters and end time from run.pipe'''

	import Compare
	import numpy as np

	stopC,stopB,stopT = False,False,False
	for i in range(nlines):
		j=nlines-i-1
#		print(A[j].split()[0])
		if ( len(A[j].split()) > 0 ):
			if (A[j].split()[0] == 'stop'):
#				print(A[j].split())
				if (time=='default'):
					indT=j
					stopT=True
				elif (A[j].split()[5] == time):
					indT=j
					stopT=True
			if (stopT == True):
				if (  A[j].split()[0] == 'aC'):
#					print(A[j].split())
					indC=j
					stopC=True
				elif (A[j].split()[0] == 'aB'):
#					print(A[j].split())
					indB=j
					stopB=True
			if (stopC & stopB & stopT):
				break
			if (i==nlines-1):
				print('Incomplete matching, stops = {0} {1} {2}'.format(
						stopC,stopB,stopT))
	if (stopC & stopB & stopT):
		aeiB = np.array(A[indB].split())[ [2,5,8] ]
		aeiC = np.array(A[indC].split())[ [2,5,8] ]
		T = A[indT].split()[ 5 ]
	else:
		stopC,stopB,stopT,indT,indB,indC = Compare.BackupParams(
									A,nlines,time)
		if np.isnan(indB):
			aeiB = ['nan', 'nan','nan']
		else:
                	aeiB = np.array(A[indB].split())[ [2,5,8] ]
		if np.isnan(indC):
			aeiC = ['nan', 'nan','nan']
		else:
                	aeiC = np.array(A[indC].split())[ [2,5,8] ]
                T = float(A[indT].split()[4])/365.25

	return (aeiB,aeiC,T)
###############################################################################
def BackupParams(A,nlines,time='default'):
	'''Backup method to find orbital parameters for B-2 sims'''

	import numpy as np

	stopT,stopB,stopC = False,False,False
	indT, indB, indC  = 0,0,0
	WPind = [nlines-1]
	for i in range(nlines):
		j=nlines-i-1
		if ( len(A[j].split()) > 0 ):
			if ((A[j].split()[0] == 'WriteParam,') | 
			    (A[j].split()[0] == 'writeparam.bash,')):
				WPind.insert(0,j)
				t = float(A[j].split()[4])/365.25
				if (time=='default'):
					print('found correct WriteParam')
					break
				elif (float(time)==t):
					print('found correct WriteParam')
					break
	print(WPind)
	indT = WPind[0]
	for i in range(WPind[1],WPind[0]-1,-1):
		if (A[i].split()[0] == 'aB'):
			if (stopB==False):
				indB=i
				stopB=True
		if (A[i].split()[0] == 'aC'):
			if (stopC==False):
				indC=i
				stopC=True
		if (stopB & stopC):
			break
		if (i==WPind[0]):
			print('end reached, B or C params missing')
			if (stopB==False):
				indB=float('NaN')
			if (stopC==False):
				indC=float('NaN')
	return(stopC,stopB,stopT,indT,indB,indC)
###############################################################################
