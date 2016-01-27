
import numpy as np
import AlphaCenModule as AC
import matplotlib.pyplot as plt
from astropy.io import ascii


### Laptop data locations
Dirs1 = ['d'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['a'+   '{0}'.format(i) for i in range(1,13)]
Dirs3 = ['b'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['c'+'{0:02}'.format(i) for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4
#Dirs = Dirs[0:2]
ndir = len(Dirs)
minpC = np.array([ ])
summary=np.zeros([0,5],dtype=[('loc', 'S3'),('vers', 'i1'), ('minpC', 'f10'), ('t', 'f10'), ('iMf', 'f10')])

for i,d in enumerate(Dirs):
### Read in TimeData.txt
	files = ['TimeData/TimeData-'+d+'-O.txt',
			 'TimeData/TimeData-'+d+'-2.txt',
			 'TimeData/TimeData-'+d+'-3.txt']
	data = [np.loadtxt(f, skiprows=1) for f in files]
	header = [open(f ,'r').readline().split() for f in files]
	for l in range(len(data)):
		assert len(header[l]) == data[l].shape[1],'Header length and # columns do not match!'
### Put this sim's data into a dictionary
	dct = [{},{},{}]
	for l in range(len(data)):
		for j,h in enumerate(header[l]):
			dct[l][h] = data[l][:,j]
		dct[l]['apB'] = dct[l]['aB']*(1.+dct[l]['eB'])
		dct[l]['apC'] = dct[l]['aC']*(1.+dct[l]['eC'])
		### Get C's closest approach, save to array

	t   = [dct[l]['t']   for l in range(len(dct))]
	aC  = [dct[l]['aC']  for l in range(len(dct))]
	eC  = [dct[l]['eC']  for l in range(len(dct))]
	iC  = [dct[l]['iC']  for l in range(len(dct))]
	pC  = [dct[l]['pC']  for l in range(len(dct))]
	apC = [dct[l]['apC'] for l in range(len(dct))]
	aB  = [dct[l]['aB']  for l in range(len(dct))]
	eB  = [dct[l]['eB']  for l in range(len(dct))]
	iB  = [dct[l]['iB']  for l in range(len(dct))]
	pB  = [dct[l]['pB']  for l in range(len(dct))]
	apB = [dct[l]['apB'] for l in range(len(dct))]
	iM  = [ np.abs(iC[l] - iB[l]) for l in range(len(iC))]
#	mgas = 10. * np.array([ max(1.e-2, (4.4e5-ti)/4.4e5) for ti in t ])

##### Make plot
	f, ax = plt.subplots(4, sharex=True)
	axtwin = ['','','']
	for ind,l in enumerate([0,2]):
		ax[ind].set_xscale('log')
		ax[ind].set_yscale('log')
		ax[ind].text(0.9, 0.95, 'a       e',
	        verticalalignment='top', horizontalalignment='right',
	        transform=ax[ind].transAxes,
	        color='black', fontsize=10)
	### Outer binary
		if not any(np.isnan(aC[l])):
			ax[ind].scatter(t[l][ aC[l] > 0.], aC[l][ aC[l] > 0.],c='black',lw=0)
			ax[ind].fill_between(t[l], pC[l], apC[l], where=apC[l] > 0.,  
												facecolor='blue', alpha=0.4)
			(dummy,y2) = ax[ind].get_ylim()
			ax[ind].fill_between(t[l], pC[l], y2, where=apC[l] <= 0.,  
												facecolor='grey', alpha=0.4)
			ax[ind].text(0.9, 0.85, '{0: 6.0f}  {1:.2f}'.format(aC[l][-1],eC[l][-1]),
		        verticalalignment='top', horizontalalignment='right',
		        transform=ax[ind].transAxes,
		        color='blue', fontsize=10)
			# calculate minimum pericenter, last time it happens, and final i
			minpC  = np.min(pC[l])
			minloc = np.where(pC[l] == minpC)
			tmin   = t[l][minloc]
			iMf    = abs( dct[l]['iC'][-1]-dct[l]['iB'][-1] )
			newline = np.array([(d, l, minpC, max(tmin), iMf)],dtype=summary.dtype)
			summary = np.append(summary,newline)
			ax[ind].plot(tmin,[minpC for ti in tmin], 'oy')
			ax[ind].text(0.05, 0.9, 'min(pC) = {0:.1f}'.format(minpC),
		        verticalalignment='top', horizontalalignment='left',
		        transform=ax[ind].transAxes,
		        color='blue', fontsize=10)
			ax[ind].text(0.05, 0.8, 'final iM = {0:.2f}'.format(iMf),
		        verticalalignment='top', horizontalalignment='left',
		        transform=ax[ind].transAxes,
		        color='blue', fontsize=10)
	### Inner binary
		ax[ind].scatter(t[l],aB[l],c='black',lw=0)
		ax[ind].fill_between(t[l], pB[l], apB[l], facecolor='red', alpha=0.4)
		ax[ind].text(0.9, 0.75, '{0: 6.1f}  {1:.2f}'.format(aB[l][-1],eB[l][-1]),
	        verticalalignment='top', horizontalalignment='right',
	        transform=ax[ind].transAxes,
	        color='red', fontsize=10)
	### Finish plot
		(y1,dummy) = ax[ind].get_ylim()
		if l == 1:
			y2 = dummy
		ax[ind].set_ylim((y1,y2))
		ax[ind].vlines(x=440000.,ymin=y1,ymax=y2)
		ax[ind].set_xlabel('Time (yrs)')
		ax[ind].set_ylabel('Distance (AU)')
	### Plot i over time
	ax[2].set_xscale('log')
	ax[2].plot(t[0], iB[0], 'r-')
	if not any(np.isnan(aC[0])):
		ax[2].plot(t[0], iM[0], 'k-',
					 t[0], iC[0], 'b-')
	ax[2].plot(t[2], iB[2], 'r--')
	if not any(np.isnan(aC[2])):
		ax[2].plot(t[2], iM[2], 'k--',
					 t[2], iC[2], 'b--')
	ax[2].set_xlabel('t (yr)')
	ax[2].set_ylabel('i (deg)')
	# Plot e
	ax[3].set_xscale('log')
	ax[3].plot(t[0], eB[0], 'r-')
	if not any(np.isnan(aC[0])):
		ax[3].plot(t[0], eC[0], 'b-')
	ax[3].plot(t[2], eB[2], 'r--')
	if not any(np.isnan(aC[2])):
		ax[3].plot(t[2], eC[2], 'b--')
	ax[3].set_ylabel('Eccentricity')

	ax[1].yaxis.tick_right()
	ax[1].yaxis.set_ticks_position('both')
	ax[1].yaxis.set_label_position("right")
	ax[3].yaxis.tick_right()
	ax[3].yaxis.set_ticks_position('both')
	ax[3].yaxis.set_label_position("right")

	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	plt.savefig('TimeData/AvsT-'+d+'.png')
	plt.clf()
	plt.close(f)

ascii.write(summary,'TimeData/PeriSummary.txt',
	names=['loc','vers','minpC','  t','  iMf'])

minp  = summary['minpC']
minpO = minp[summary['vers']==0]
minp3 = minp[summary['vers']==2]

stats = np.array( [ [np.median(l) for l in [minp, minpO, minp3]],
					[np.mean(l)   for l in [minp, minpO, minp3]],
					[np.min(l)    for l in [minp, minpO, minp3]],
					[np.max(l)    for l in [minp, minpO, minp3]] ])
print(stats)
### Plot pericenters
f2, ax2 = plt.subplots(2,2)

ax2[0,0].set_xscale('log')
#ax2[l].set_yscale('log')
ax2[0,0].scatter( summary['minpC'],summary['iMf'] )
ax2[0,1].hist(minp)
ax2[0,1].set_xlabel('All min(pericenter)')
ax2[1,0].hist(minpO)
ax2[1,0].set_xlabel('Orig min(pericenter)')
ax2[1,1].hist(minp3)
ax2[1,1].set_xlabel('Trip min(pericenter)')


plt.savefig('TimeData/Summaries.png')
plt.clf()
plt.close(f2)




