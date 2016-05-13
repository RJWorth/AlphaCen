### Read TimeData.txt files from all simulations and compile into DiskSummary.txt
### Pipeline: simulations run on hammer
#             important files copied onto chloe's local director
#             WriteTimeData.py runs AlphaCenModule.WriteTimeData on each directory
#             Compare.bash copies the TimeData.txt files from chloe to laptop
#             ReadAllTimeDatas.py reads the TimeData.txt files and creates DiskSummary.txt
#             ReadDiskSum.py reads DiskSummary.txt and does analysis on it

import numpy as np
import AlphaCenModule as AC
import Merc
import matplotlib.pyplot as plt
from astropy.io import ascii
from statsmodels.robust.scale import mad as madfn
import pickle
import pandas as pd
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

### Laptop data locations
Dirs1 = ['d'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['a'+   '{0}'.format(i) for i in range(1,13)]
Dirs3 = ['b'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['c'+'{0:02}'.format(i) for i in range(1,9)]

#Dirs = Dirs2+Dirs3+Dirs4+Dirs1
Dirs = Dirs4+Dirs1+Dirs2+Dirs3
Dirs = ['b24','d15']
ndir = len(Dirs)

### Code for Original, DiskB-2, and DiskB-3 versions of each sim
subdirs = ['O','A2','A3','B2','B3']

### Create empty array to store rows of data in as I iterate over the directories
cols =  ['dr','vrs',  'tf',  
           'aBf', 'eBf', 'pBf',  'apBf',
           'aCf', 'eCf', 'pCf',  'apCf',
           'minpB','mintB','minpC','mintC', 
           'iMf','rtrf','rpf','rpm']

#-----------------------------------------------------------------------------#

for i, d in enumerate(Dirs):
	for j, s in enumerate(subdirs):
		# A dirs don't exist for sets a and b
		if not (any([letter in d for letter in ['a','b']]) & ('A' in s)):
			simID = d+'-'+s
			thisdir = pd.read_csv('TimeData/TimeData-'+simID+'.txt',delim_whitespace=True)
			# Calculate non-time-dependent stats for this sim, and add to summary
			dr    = d
			vrs   = s
			tf    = thisdir.t.iloc[-1]
			aBf   = thisdir.aB.iloc[-1]
			eBf   = thisdir.eB.iloc[-1]
			pBf   = aBf*(1-eBf)
			apBf  = aBf*(1+eBf)
			aCf   = thisdir.aC.iloc[-1]
			eCf   = thisdir.eC.iloc[-1]
			pCf   = aCf*(1-eCf)
			apCf  = aCf*(1+eCf)
			minpB = min(thisdir.pB)
			mintB = max(thisdir.t.iloc[Merc.which(thisdir.pB,minpB)])
			if not all(np.isnan(thisdir.pC)):
				minpC = min(thisdir.pC)
				mintC = max(thisdir.t.iloc[Merc.which(thisdir.pC,minpC)])
			elif '2' in s:
				minpC, mintC = np.nan, np.nan
			else:
				print('Warning: pC = nan, but subdir is {}, not a disk binary!'.format(s))
			iMf   = thisdir.iC.iloc[-1]-thisdir.iB.iloc[-1]
			rtrf  = thisdir.rtr.iloc[-1]
			rpf   = rtrf/pBf
			rpm   = rtrf/minpB

			# make a dictionary of these stats, to enter as one line in the dataframe
			line = {}
			for c in cols:
				line[c] = eval(c)

			### Add line to dataframe
			# if first time through, create summary frame
			if (i==0) and (j==0):
				counter = 0
				summary = pd.DataFrame( line, index=[0])
				summary = summary[cols]
			# otherwise add to the dataframe
			else:
				counter += 1
				summary.loc[counter] = line

with open('TimeData/DiskSummary.txt','w') as f:
	f.write(summary.to_string(justify='right'))



