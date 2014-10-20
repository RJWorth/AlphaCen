
def GetSubdirs(indir):
	import os
	return [name for name in os.listdir(indir)
			if os.path.isdir(os.path.join(indir, name))]


def NextSubDir(searchdir,prefix):
	import numpy as np
	from math import ceil

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

### Report outcome
	print('Dirs:   {}'.format(dirs))
	print('Use:    {}'.format(ValidInd))
	print('Ignore: {}'.format(InvalidInd))

	print('Max directory number found:  {}'.format(MaxUsed))
	print("Numerals' median digits:     {}".format(ndigits))

	print('First open directory number: {}'.format(OpenDir))

	Name = prefix + str(OpenDir).zfill(ndigits)
	print('Resulting name: {}'.format(Name))
	return(Name)

