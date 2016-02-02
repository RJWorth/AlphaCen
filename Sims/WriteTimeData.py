
import numpy as np
import AlphaCenModule as AC

### Chloe/hammer data locations
Dirs1 = [       'Proxlike/Prx'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['Proxlike/071714/Prx'+  '{0:}'.format(i) for i in range(1,13)]
Dirs3 = ['Proxlike/072314/Prx'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['Proxlike/081414/Prx'+'{0:02}'.format(i) for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4

#Dirs = Dirs[0:2]

for i,d in enumerate(Dirs):
	thisdir = d+'/Original'
	print(thisdir)
	AC.WriteAEI(thisdir)	

	thisdir = d+'/DiskB-2'
	print(thisdir)
	AC.WriteAEI(thisdir)	

	thisdir = d+'/DiskB-3'
	print(thisdir)
	AC.WriteAEI(thisdir)	




