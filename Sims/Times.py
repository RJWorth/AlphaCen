def GetTime(b):
	'''Simulation time derived from bytes in xv.out'''

### Modules
	import numpy as np

### Fit parameters from earlier R calculation
	yIntercept = -126072.
	slope      =   66773.

	t = 10.**( (b - yIntercept)/slope )

	return(t)

###############################################################################
def WallTime(simt, machine):
	'''Wall time derived from simulation time, based on models'''

### Setup
	assert (machine == 'chloe') | (machine == 'shapiro')

### Modules
	import numpy as np

	if   (isinstance( simt,np.ndarray )):
		flag='np'
	elif (isinstance( simt,      list )):
		flag='arr'
	elif (isinstance( simt,     float )):
		flag='flt'
	else:
		print('From WallTime: What type is simt?')
	
	if   flag=='flt':
		simt=np.array([simt])
	elif flag=='arr':
		simt=np.array( simt )

### Read coefficients
	cfile = open('RunTimeCoefs.txt')
	coefs = cfile.readlines()
	cfile.close()

	coefs = np.array([row.split()[1:] for row in coefs[1:]])
	coefs = coefs.astype(float)

### Loop over items in simt
	WallT = np.empty_like(simt)
	for i,t in enumerate(simt):
### Determine which coefficients to use
		if   (machine ==   'chloe') & (t  < 1e6):
			row = 0
		elif (machine ==   'chloe') & (t >= 1e6):
			row = 2
		if   (machine == 'shapiro') & (t  < 1e6):
			row = 1
		elif (machine == 'shapiro') & (t >= 1e6):
			row = 3
	
		c, exp = coefs[row,:]

### Fit parameters from earlier R calculation
		if t < 1e6:
			WallT[i] = c**(t**exp)	# double exponential model
		else:
			WallT[i] = c *(t**exp)	# regular exponential model

	if flag=='flt':
		WallT=WallT[0]

	return(WallT)

###############################################################################
def PredictBytes(t):
	'''Bytes as a function of time'''

### Modules
	import numpy as np
	from math import log10

### Fit parameters from earlier R calculation
	yIntercept = -126072.
	slope      =   66773.

	bytes = yIntercept + slope * log10(t)

	return(bytes)

###############################################################################
def TimeRemaining(Dir, machine='default'):

### Modules
	import os 
	import Times
	import datetime as DT
	from math import log10, floor

### Constants
	MaxSize = 474885			# bytes in 1e9-yr simulation's xv.out file
### b[i] = bytes in 1e[i] sim
	b = [Times.PredictBytes(10.**i) for i in range(0,10)]

### Variables
	if (machine == 'default'):
		if (Dir[0] == 'S'):
			machine='shapiro'
		elif (Dir[0] == 'C'):
			machine='chloe'
		else:
			machine='chloe'
			print('Warning: unrecognized machine, using chloe models')

	assert (machine=='shapiro') | (machine=='chloe')

	inTime = os.path.getmtime(Dir+'/In/big.in')		# start of whole run
	rpTime = os.path.getmtime(Dir+'/run.pipe')		# start of current step
	xvTime = os.path.getmtime(Dir+'/Out/xv.out')	# now/most recent write

	xvSize = os.path.getsize( Dir+'/Out/xv.out')	# number of bytes

### Calculations
	# sim time now, and log(simtime) at start and end of this timestep
	CurrentYear    = Times.GetTime(xvSize)
	Last		   = int(floor(log10(CurrentYear)))
	Next		   = int(floor(log10(CurrentYear))+1.)

	# measured elapsed times
	ElapsedTime_a = xvTime-inTime				# time elapsed (all)
	ElapsedTime_s = xvTime-rpTime				# time elapsed (this step)

	# initial and final times for this timestep, and current time
	ExpectedWallTi = Times.WallTime(   10.**Last, machine)
	ExpectedWallTc = Times.WallTime( CurrentYear, machine)
	ExpectedWallTf = Times.WallTime(   10.**Next, machine)
	ExpectedWallT9 = Times.WallTime(     10.**9., machine)

	# corresponding actual times
	ElapsedWallTi  = rpTime-inTime		# walltime to start of step
	ElapsedWallTc  = xvTime-inTime		# walltime to current time
	ElapsedWallTs  = xvTime-rpTime		# walltime elapsed this step only

	# Byte calculations
	CompletedFraction = float(xvSize-b[Last])/float(b[Next]-b[Last])
	BytesRemaining = b[Next]-xvSize
	BytesRemaining9= b[9]   -xvSize
	ByteProgress   = xvSize-b[Last]				# subtract size at 1e8 yrs
	Rate 		   = ByteProgress/ElapsedTime_s	# avg bytes/second
	TimeRemaining  = BytesRemaining /Rate	
	TimeRemaining9 = BytesRemaining9/Rate	

### Get expected time of completion with old method
	# for whole 1e9 sim
	t   = DT.datetime.fromtimestamp(xvTime)
	dt9 = DT.timedelta(seconds=TimeRemaining9)
	tf9 = t+dt9
	# for this timestep
	if (Next < 9):
		dt  = DT.timedelta(seconds=TimeRemaining)
		tf  = t+dt

	ExpectedTime_a = ExpectedWallTc
	if (ElapsedTime_a > ExpectedTime_a):
		comparison = 'behind'
	else:
		comparison = 'ahead'

### Get expected time of completion with new method
	# for whole 1e9 sim
	dt9_new = DT.timedelta(seconds=ExpectedWallT9-ExpectedWallTc)
	tf9_new = t+dt9_new
	# for this timestep	
	if (Next < 9):
		dt_new  = DT.timedelta(seconds=ExpectedWallTf-ExpectedWallTc)
		tf_new  = t+dt_new

### Write information to terminal
#	print('{0}: {1}/{2} bytes, {3:.1f}% complete, t = {4:.3g}.'.format(
#			Dir,xvSize,MaxSize, CompletedFraction*100.,CurrentYear))
	print('{0:>6}: {1:8.3g} {2:.1f}%'.format(
			Dir, CurrentYear, CompletedFraction*100.)+
		  '  {0:6.0f}/{1:6.0f} {2:>6}'.format(
			ElapsedTime_a, ExpectedTime_a, comparison)+
		  '    {0} / {1}'.format(
				  str(tf9)[5:16], str(tf9_new)[5:16] ))
	if (Next < 9):
		print('                                           1e{0}: {1} / {2}'.format(
			Next, str(tf)[5:16],  str(tf_new)[5:16] ))

#	print('    Next:  Min  /   Hr  /    D    to completion')
#	print('    1e{0} {1:6.1f}  {2:6.2f}  {3:6.2f}'.format(
#		Next,TimeRemaining/60.,TimeRemaining/3600.,TimeRemaining/(24.*3600.)))

