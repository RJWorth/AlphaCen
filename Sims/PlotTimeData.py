
import numpy as np
import AlphaCenModule as AC
import matplotlib.pyplot as plt

### Chloe data locations
Dirs1 = [       'Proxlike/Prx'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['Proxlike/071714/Prx'+  '{0:}'.format(i) for i in range(1,13)]
Dirs3 = ['Proxlike/072314/Prx'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['Proxlike/081414/Prx'+'{0:02}'.format(i) for i in range(1,9)]

### Laptop data locations
Paths = Dirs1+Dirs2+Dirs3+Dirs4
Paths = Dirs1+Dirs2+Dirs4 # no subdirs in Dirs3

for i,d in enumerate(Dirs):
	thisdir = d+'/DiskB-2'
	print(thisdir)
	AC.WriteAEI(thisdir)	


