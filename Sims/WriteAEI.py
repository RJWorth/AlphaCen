
import numpy as np
import AlphaCenModule as AC

Dirs1 = [       'Proxlike/Prx'+'{0:02}'.format(i)+'/Original' for i in range(1,35)]
Dirs2 = ['Proxlike/071714/Prx'+  '{0:}'.format(i)+'/Original' for i in range(1,13)]
Dirs3 = ['Proxlike/072314/Prx'+'{0:02}'.format(i)             for i in range(1,27)]
Dirs4 = ['Proxlike/081414/Prx'+'{0:02}'.format(i)+'/Original' for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4
ndir = len(Dirs)

print(Dirs)
print(Dirs[0])
for d in [Dirs[0]]:
### Create TimeData.txt file for each, in Out/Aei...
##	AC.WriteAEI(d)

### Read in TimeData.txt
	timedata = np.loadtxt(d+'/Out/AeiOutFiles/TimeData.txt', skiprows=1)
	header = open(d+'/Out/AeiOutFiles/TimeData.txt' ,'r').readline().split()

	print(timedata)
	print(header)
	print(timedata.shape)

### Get C's closest approach, save to array
	


