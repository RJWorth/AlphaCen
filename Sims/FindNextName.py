
def GetSubdirs(indir):
	import os
	return [name for name in os.listdir(indir)
			if os.path.isdir(os.path.join(indir, name))]

###############################################################################
def NextSubDir(searchdir,prefix='default',PrintInfo=False):
	import numpy as np
	from math import ceil

### Setup prefix
	if prefix=='default':
		if searchdir=='Proxlike':
			prefix='Prx'
		elif searchdir=='Err':
			prefix=''
		else:
			print('Prefix unspecified with no relevant default')
			assert 1==0

	dirs = np.array(GetSubdirs(searchdir))

### Index for valid dirs, and array of valid dirs
	ValidInd   = np.array([ (prefix in d) & 
							d.replace(prefix,'').isdigit() 
								for d in dirs])
	InvalidInd = np.array([not i for i in ValidInd])

	ValidStr   = np.array(dirs[ValidInd])
	InvalidStr = np.array(dirs[InvalidInd])
### Convert valid dir numerals to integers
	ValidNumerals = np.array( [d.replace(prefix,'') for d in ValidStr] )
	ValidNumeric  = np.array( [int(d) for d in ValidNumerals] )

	MaxUsed = max(ValidNumeric)

### Determine correct number of digits
	ndigits = int(ceil(np.median( [len(d) for d in ValidNumerals] )))

### Iterate from 1 up to the maximum already-used integer, until reaching empty slot
	for i in range(1,MaxUsed+2):
		if i not in ValidNumeric:
			OpenDir = i
			break

	Name = prefix + str(OpenDir).zfill(ndigits)

### Report outcome
	if PrintInfo == True:
		print('Dirs:   {0}'.format(dirs))
		print('Use:    {0}'.format(ValidInd))
		print('Ignore: {0}'.format(InvalidInd))
		print('Max directory number found:  {0}'.format(MaxUsed))
		print("Numerals' median digits:     {0}".format(ndigits))
		print('First open directory number: {0}'.format(OpenDir))
		print('Resulting name: {0}'.format(Name))

	return(Name)

