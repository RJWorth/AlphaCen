###############################################################################
def DiskFateSum(WhichDirs='default'):
	'''Get percentage of disk material that met each fate'''

	import AlphaCenModule as AC
	import numpy as np

	if WhichDirs=='default':
		WhichDirs = ['Proxlike/Prx{0:02d}'.format(i) for i in range(1,35)]
	
	for dr in WhichDirs:
		a, dest2, b, c, d = AC.ReadInfo(dr+'/DiskB-2')
		a, dest3, b, c, d = AC.ReadInfo(dr+'/DiskB-3')
		dest2,dest3=np.array(dest2),np.array(dest3)
		fates2 = np.array([len(dest2),sum(dest2=='ejected'),
				sum(dest2=='AlCenB'),sum(dest2=='Center'),sum(dest2=='PrxCen')])
		fates3 = np.array([len(dest3),sum(dest3=='ejected'),
				sum(dest3=='AlCenB'),sum(dest3=='Center'),sum(dest3=='PrxCen')])
		print(dr+'    {0} {1}     {2} {3}'.format(sum(fates2[1:]),fates2,
							  sum(fates3[1:]),fates3 ))


