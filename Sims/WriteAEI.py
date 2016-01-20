
import numpy as np
import AlphaCenModule as AC
import matplotlib.pyplot as plt


Dirs1 = [       'Proxlike/Prx'+'{0:02}'.format(i)+'/Original' for i in range(1,35)]
Dirs2 = ['Proxlike/071714/Prx'+  '{0:}'.format(i)+'/Original' for i in range(1,13)]
Dirs3 = ['Proxlike/072314/Prx'+'{0:02}'.format(i)             for i in range(1,27)]
Dirs4 = ['Proxlike/081414/Prx'+'{0:02}'.format(i)+'/Original' for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4
ndir = len(Dirs)
minpC = np.array(len(Dirs))

for i,d in enumerate([Dirs[0]]):
### Create TimeData.txt file for each, in Out/Aei...
##	AC.WriteAEI(d)

### Read in TimeData.txt
	timedata = np.loadtxt(d+'/Out/AeiOutFiles/TimeData.txt', skiprows=1)
	header = open(d+'/Out/AeiOutFiles/TimeData.txt' ,'r').readline().split()

	print(header)
	assert len(header) == timedata.shape[1],'Header length and # columns do not match!'

	dict = {}
	for j,h in enumerate(header):
		dict[h] = timedata[:,j]

### Get C's closest approach, save to array
	minpC[i] = np.min(dict['pC'])	
	print(dict.items())

	minloc = np.where(dict['pC'] == minpC[i])
	print(dict['pC'][minloc])
	print(dict['t'][minloc])
	print(dict['pB'][minloc])

	##### Make plot
	f, ax1 = plt.subplots(1)
	ax.scatter(dict['t'],dict['pC'])
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlabel('Time')
	ax1.set_ylabel('Pericenter')
	  	
	plt.savefig(d+'/PvsT.pdf')
	plt.clf()
	plt.close(f)



