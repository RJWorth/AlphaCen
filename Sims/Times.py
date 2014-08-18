def GetTime(bytes):
	'''Simulation time derived from bytes in xv.out'''

### Fit parameters from earlier R calculation
	yIntercept = -126072.
	slope      = 66773.

	t = 10.**( (bytes - yIntercept)/slope )

	return(t)

###############################################################################
def PredictBytes(t):
	'''Bytes as a function of time'''

	from math import log10

### Fit parameters from earlier R calculation
	yIntercept = -126072.
	slope      = 66773.

	bytes = yIntercept + slope * log10(t)

	return(bytes)

###############################################################################
def TimeRemaining(Dir):

	import os 
	import Times
	import datetime as DT
	from math import log10, floor

### Constants
	MaxSize = 474885			# bytes in 1e9-yr simulation's xv.out file
### b[i] = bytes in 1e[i] sim
	b = [Times.PredictBytes(10.**i) for i in range(0,10)]

### Variables
	rpTime = os.path.getmtime(Dir+'/run.pipe')
	xvTime = os.path.getmtime(Dir+'/Out/xv.out')

	xvSize = os.path.getsize( Dir+'/Out/xv.out')	# number of bytes

### Calculations
	CurrentYear    = Times.GetTime(xvSize)
	Last		   = int(floor(log10(CurrentYear)))
	Next		   = int(floor(log10(CurrentYear))+1.)

	ElapsedTime = xvTime-rpTime				# time elapsed from 1e8 yrs to now
	CompletedFraction = float(xvSize-b[Last])/float(b[Next]-b[Last])


	BytesRemaining = b[Next]-xvSize
	ByteProgress   = xvSize-b[Last]				# subtract size at 1e8 yrs
	Rate 		   = ByteProgress/ElapsedTime	# avg bytes/second
	TimeRemaining  = BytesRemaining/Rate	

	t  = DT.datetime.fromtimestamp(xvTime)
	dt = DT.timedelta(seconds=TimeRemaining)
	tf = t+dt

### Write information to terminal
	print('{0}: {1}/{2} bytes, {3:.1f}% complete, t = {4:.3g}.'.format(
			Dir,xvSize,MaxSize, CompletedFraction*100.,CurrentYear))
	print('    Average rate: {0:.3f} bytes/second.'.format(Rate))
	print('    Expected time of completion: {0}'.format(tf))
	print('    Next:  Min  /   Hr  /    D    to completion')
	print('    1e{0} {1:6.1f}  {2:6.2f}  {3:6.2f}'.format(
		Next,TimeRemaining/60.,TimeRemaining/3600.,TimeRemaining/(24.*3600.)))

