###############################################################################
def ComparePipedParams(WhichDir):
	'''Compare orbital parameters'''
	
### Modules
	import AlphaCenModule as AC
	import numpy as np
	import Compare

### Read in the three relevant run.pipe files
	nOlines = AC.FileLength(WhichDir+'/Original/run.pipe')
	origfile=open(          WhichDir+'/Original/run.pipe')
	O=origfile.readlines()
	origfile.close()

	nAlines = AC.FileLength(WhichDir+'/DiskA-3/run.pipe')
	Afile=open(          WhichDir+'/DiskA-3/run.pipe')
	A=Afile.readlines()
	Afile.close()

	nBlines = AC.FileLength(WhichDir+'/DiskB-3/run.pipe')
	Bfile=open(          WhichDir+'/DiskB-3/run.pipe')
	B=Bfile.readlines()
	Bfile.close()

### Check lengths against each other
	if (nAlines != nBlines):
		print('Warning: lengths of run.pipe for A-3 and B-3 do not match!')

### Find param lines from end of A
	AaeiB,AaeiC,AT = Compare.FindParams(A,nAlines)
	BaeiB,BaeiC,BT = Compare.FindParams(B,nBlines)
	if (AT != BT):
		print('Warning: lengths of run.pipe for A-3 and B-3 do not match! {0} {1}'.format(AT,BT))
	if ( max(float(AT),float(BT)) == float(AT) ):
		T = AT
	else:
		T = BT
	OaeiB,OaeiC,OT = Compare.FindParams(O,nOlines,T)

	Binary = np.array([OaeiB,AaeiB,BaeiB])
	Triple = np.array([OaeiC,AaeiC,BaeiC])
	allT = np.array([OT,AT,BT])

	display=np.c_[['O','A','B'],Binary,Triple,allT]
	print(['','   aB','    eB','   iB','   aC','   eC','   iC','t'])
	print(display)

	names=['a','e','i']
	if ((allT[0]!=allT[1]) | (allT[0]!=allT[2])):
		print('Mismatch in time: {0}'.format(allT))
	for i in range(3):
		if ((Binary[0,i]!=Binary[1,i]) | (Binary[0,i]!=Binary[2,i])):
			print('***Mismatch in {0}: {1}***'.format(names[i],Binary[:,i]))


###############################################################################
def FindParams(A,nlines,time='default'):
	'''Get B and C orbital parameters and end time from run.pipe'''

	import numpy as np

	stop1,stop2,stop3 = False,False,False
	for i in range(nlines):
		j=nlines-i-1
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
	aeiB = np.array(A[indB].split())[ [2,5,8] ]
	aeiC = np.array(A[indC].split())[ [2,5,8] ]
	T = A[indT].split()[ 5 ]

	return (aeiB,aeiC,T)
