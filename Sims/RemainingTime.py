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
	import RemainingTime as RT

	rpTime = os.path.getmtime(Dir+'/run.pipe')
	xvTime = os.path.getmtime(Dir+'/Out/xv.out')

	xvSize = os.path.getsize( Dir+'/Out/xv.out')
	MaxSize = 474885			# bytes in 1e9-yr simulation's xv.out file
	b8 = RT.PredictBytes(1e8)	# bytes in 1e8 sim
	b9 = RT.PredictBytes(1e9)	# bytes in 1e9 sim
	assert b9==MaxSize


	ElapsedTime = xvTime-rpTime		# time elapsed from 1e8 yrs to now

# fraction of xv.out size written
	CompletedFraction = float(xvSize-b8)/float(b9-b8)
	CurrentYear = RT.GetTime(xvSize)

	BytesRemaining = b9-xvSize
	ByteProgress   = xvSize-b8					# subtract size at 1e8 yrs
	Rate 		   = ByteProgress/ElapsedTime	# avg bytes/second
	TimeRemaining  = BytesRemaining/Rate	

	print('{0}: {1}/{2} bytes, {3:.1f}% complete, year {4:.3g}.'.format(
			Dir,xvSize,MaxSize, CompletedFraction*100.,CurrentYear))
	print('Average rate: {0:.3f} bytes/second.'.format(Rate))
	print('Estimated time remaining: {0:.1f} minutes.'.format(
		TimeRemaining/60.))
	print('                          {0:.2f} hours.  '.format(
		TimeRemaining/3600.))
	print('                          {0:.2f} days.   '.format(
		TimeRemaining/(24.*3600.)))

