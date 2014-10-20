
import os
import numpy as np
import AlphaCenModule as AC

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

### Get names of err directories and which ones have solutions
ErrDirs = np.array(sorted(get_immediate_subdirectories('Err')))
solved  = np.array([len(d)>2 for d in ErrDirs])

f = open('Err/'+dirs[0]+'/summary.out','r')
summ = f.readlines()
lastlines = [np.array(np.array(summ[-1].split()))]
for d in dirs[solved][1:sum(solved)]:
	f = open('Err/'+d+'/summary.out','r')
	summ = f.readlines()
	lastlines = np.append(lastlines,
						[np.array(summ[-1].split())],axis=0)

### Names of regular simulation directories
SimDirs = ['SDir'+str(i) for i in range(1,8)] + \
		  ['CDir'+str(i) for i in range(1,11)]

summatch = np.array([[i,'',np.nan] for i in ErrDirs[solved]])
for i,d in enumerate(SimDirs):
	f = open(d+'/summary.out','r')
	s = f.readlines()
	s = np.array([line.split() for line in s])
	for i,row1 in enumerate(lastlines):
		for j,row2 in enumerate(s):
			if sum(row1==row2)>10:
				print('Match {0}: lastlines[{1}] = {2}[{3}]'.format(
														sum(row1==row2),i,d,j))
				summatch[i,1:3]=[d,j]

summatch[summatch[:,1].argsort()]	# in order of SimDir


### Write new summary line to fixed errdirs -- DO ONLY ONCE
for d in ErrDirs[solved][1:]:
	AC.Summary('Err/'+d,1e9,WhichTime='Err',mA=.123,mB=.123)



