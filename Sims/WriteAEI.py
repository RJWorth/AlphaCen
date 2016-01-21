
import numpy as np
import AlphaCenModule as AC
import matplotlib.pyplot as plt

### Laptop data locations
Dirs1 = ['d'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['a'+'{0:02}'.format(i) for i in range(1,13)]
Dirs3 = ['b'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['c'+'{0:02}'.format(i) for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4
#Dirs = [Dirs[0]]
ndir = len(Dirs)
minpC = np.array([0. for i in range(ndir)])

for i,d in enumerate(Dirs):
### Read in TimeData.txt
	datafile = 'TimeData/TimeData-'+d+'.txt'
	timedata = np.loadtxt(datafile, skiprows=1)
	header = open(datafile ,'r').readline().split()
	print(header)
	assert len(header) == timedata.shape[1],'Header length and # columns do not match!'
### Put this sim's data into a dictionary
	ddata = {}
	for j,h in enumerate(header):
		ddata[h] = timedata[:,j]

	ddata['apB'] = ddata['aB']*(1.+ddata['eB'])
	ddata['apC'] = ddata['aC']*(1.+ddata['eC'])
### Get C's closest approach, save to array
	minpC[i] = np.min(ddata['pC'])
	minloc = np.where(ddata['pC'] == minpC[i])
	tmin = ddata['t'][minloc]
##### Make plot
	t   = ddata['t']
	aC  = ddata['aC']
	pC  = ddata['pC']
	apC = ddata['apC']
	mgas = 10. * np.array([ max(1.e-2, (4.4e5-ti)/4.4e5) for ti in t ])

	f, ax = plt.subplots(2) # sharex=True)
	ax[0].set_xscale('log')
	ax[0].set_yscale('log')
### Outer binary
	ax[0].scatter(t[ aC > 0.], aC[ aC > 0.],c='black',lw=0)
	ax[0].fill_between(t, pC, apC, where=apC > 0.,  facecolor='blue', alpha=0.4)
	(dummy,y2) = ax[0].get_ylim()
	ax[0].fill_between(t, pC, y2, where=apC <= 0.,  facecolor='grey', alpha=0.4)
### Inner binary
	ax[0].scatter(ddata['t'],ddata['aB'],c='black',lw=0)
	ax[0].fill_between(t, ddata['pB'], ddata['apB'], facecolor='red', alpha=0.4)
### Plot just inner, just outer binaries
	ax[1].set_xscale('log')
#	ax[1].set_yscale('log')
	ax[1].scatter(t, mgas, c='black',lw=0)

### Finish plot
	(y1,dummy) = ax[0].get_ylim()
	ax[0].set_ylim((y1,y2))
	ax[0].set_xlabel('Time (yrs)')
	ax[0].set_ylabel('Distance (AU)')
	plt.savefig('TimeData/AvsT-'+d+'.png')
	plt.clf()
	plt.close(f)



