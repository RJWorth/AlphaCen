### Analyze properties of all sims and look for correlations
### Pipeline: simulations run on hammer
#             important files copied onto chloe's local director
#             WriteTimeData.py runs AlphaCenModule.WriteTimeData on each directory
#             Compare.bash copies the TimeData.txt files from chloe to laptop
#             ReadAllTimeDatas.py reads the TimeData.txt files and creates DiskSummary.txt
#             ReadDiskSum.py reads DiskSummary.txt and does analysis on it

%matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

### Read in data
X = pd.read_csv('TimeData/DiskSummary.txt',delim_whitespace=True)
nsims = X.shape[0]

### Make boolean indices for which sims are stable, and which are Proxima-like
stb = (X.tf>=1e6) & pd.notnull(X.aBf) & (X.aBf>=0) & ((X.aCf>=0) | pd.isnull(X.aCf))
prx = stb & X.apCf>10000 & pd.notnull(X.apCf)
### Sims which have a disk
dsk = pd.notnull(X.rtrf)
### Sims with disks centered on A or B
cnA = pd.Series(['A' in s for s in X.vrs])
cnB = pd.Series(['B' in s for s in X.vrs])
### Binary or Triple disk sims
bn  = pd.Series(['2' in s for s in X.vrs])
tr  = pd.Series(['3' in s for s in X.vrs])

### binned time parameter
#tbins = [5,6,7,8]
#nbins = len(tbins)+1
#tbin = np.digitize(np.log10(X.tf), tbins, right=True)+1
#X['tbin'] = tbin

### Pair plots
cols = ['tf','aBf','eBf','apBf','pBf','minpB','rtrf']
# all stable sims, color = vrs
#sb.pairplot(X.loc[stb], vars=cols, 
#	hue='vrs', palette = ['black','orange','red','green','blue'])
# proxlike sims
sb.pairplot(X.loc[stb & cnB], vars=cols, 
	hue='vrs', palette = ['black','orange','red','green','blue'])
#sb.pairplot(X.loc[stb & cnB], vars=cols, 
#	hue='vrs', palette = sb.cubehelix_palette(5, start=.5, rot=-.5))

### Correlation matrix
X.loc[prx & cnB].corr()
X[cols].loc[prx & cnB].corr()

### Plot min vs fin pB
tricol = np.array(['black' for j in range(nsims)])
tricol[np.array(tr)] = 'red'
plt.scatter(X.pBf.loc[], X.minpB.loc[prx & cnB & bn], c=tricol[np.array(prx & cnB & bn)])

