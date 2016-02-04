
import numpy as np
import AlphaCenModule as AC
import Merc
import matplotlib.pyplot as plt
from astropy.io import ascii
from statsmodels.robust.scale import mad

### Laptop data locations
Dirs1 = ['d'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['a'+   '{0}'.format(i) for i in range(1,13)]
Dirs3 = ['b'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['c'+'{0:02}'.format(i) for i in range(1,9)]

Dirs = Dirs1+Dirs2+Dirs3+Dirs4
Dirs = Dirs[0:2]
ndir = len(Dirs)

### Code for Original, DiskB-2, and DiskB-3 versions of each sim
vers = ['O','2','3']

### Format string for columns wanted in summary
fmt =    '{0.dr: >3}'   +' {0.vrs: >3}' +' {0.tf:9.2e}'   +\
		' {0.aBf: 6.2f}'  +' {0.eBf: 8.4f}'  +' {0.aCf: 8.1f}'  +' {0.eCf: 8.4f}'  +\
		' {0.minpB:7.2f}'+' {0.mintB:9.2e}'+' {0.minpC:7.2f}'+' {0.mintC:9.2e}'+\
		' {0.iMf:6.1f}'  +' {0.rtrf:6.2f}'
### Object class to easily fill one line of data into above
#-----------------------------------------------------------------------------#
class Line(object):
	'''Line of data on one sim/vers, for summary.'''
	def __init__(self, dr,vrs,tf,
				aBf,eBf, aCf,eCf,
				minpB,mintB,minpC,mintC,
				iMf,rtrf):
		self.dr     = dr
		self.vrs   = vrs
		self.tf     = tf
		self.aBf    = aBf
		self.eBf    = eBf
		self.aCf    = aCf
		self.eCf    = eCf
		self.minpB  = minpB
		self.mintB  = mintB
		self.minpC  = minpC
		self.mintC  = mintC
		self.iMf    = iMf
		self.rtrf   = rtrf
#-----------------------------------------------------------------------------#
### header text for summary file
#hdr = Line('dr','vrs','tf', 'aBf', 'eBf', 'aCf', 'eCf',
#						'minpB','mintB','minpC','mintC',
#						'iMf','rtrf')
hdrcols = ['dr','vrs','tf', 'aBf', 'eBf', 'aCf', 'eCf', 'minpB','mintB','minpC','mintC','iMf','rtrf']
widths  = [   3,    3,   9,     6,     8,     8,     8,       7,      9,      7,      9,    6,     6]
hdrstr = ' '.join(['{0['+str(i)+']: >'+str(w)+'}' for i,w in enumerate(widths)])
hdr = hdrstr.format(hdrcols)
#summary=[hdr]
summary=[]

for i,d in enumerate(Dirs):
### Read in TimeData.txt
	files = ['TimeData/TimeData-'+d+'-'+l+'.txt' for l in vers]
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
	rtr = [dct[l]['rtr'] for l in range(len(dct))]
#	mgas = 10. * np.array([ max(1.e-2, (4.4e5-ti)/4.4e5) for ti in t ])

##### Make plots
	f, ax = plt.subplots(6, sharex=True)
	axtwin = ['','','']

	# calculate minimum pericenter, last time it happens, and final i
	for ind,l in enumerate([0,1,2]):
		# get basic stats of this sim
		dr    = d
		vrs   = vers[l]
		tf    = t[l][-1]
		aBf   = aB[l][-1]
		eBf   = eB[l][-1]
		aCf   = aC[l][-1]
		eCf   = eC[l][-1]
		minpB = min(pB[l])
		mintB = t[l][Merc.which(pB[l],minpB)]
		if not all(np.isnan(pC[l])):
			minpC = min(pC[l])
			mintC = t[l][Merc.which(pC[l],minpC)]
		elif l==1:
			minpC, mintC = np.nan, np.array([np.nan])
		else:
			print('Warning: pC = nan, but vers not B-2!')
		iMf   = iM[l][-1]
		rtrf  = rtr[l][-1]
		# add to summary
		line = Line(dr,vrs,tf,aBf,eBf,aCf,eCf,
						minpB,max(mintB),minpC,max(mintC),iMf,rtrf)
#		newline = np.array(
#					[(d, l, minpC[l], max(tmin[l]), iM[l][-1],rtr[l][-1])], 
#								dtype=summary.dtype)
#		summary = np.append(summary,newline)
		summary = np.append(summary, line)
	### Plot star orbits for Original and DiskB-3
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
			ax[ind].plot(mintC,[minpC for ti in mintC], 'oy')
			ax[ind].text(0.05, 0.9, 'min(pC) = {0:.1f}'.format(minpC),
		        verticalalignment='top', horizontalalignment='left',
		        transform=ax[ind].transAxes,
		        color='blue', fontsize=10)
			ax[ind].text(0.05, 0.8, 'final iM = {0:.2f}'.format(iM[l][-1]),
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

	# Plot rtr
	ax[4].set_xscale('log')
	ax[4].plot(t[1], rtr[1], 'r-')
	if not any(np.isnan(rtr[2])):
		ax[4].plot(t[2], rtr[2], 'b-')
		ax[4].set_ylabel('Truncation Radius (AU)')
	(y4a, y4b) = ax[4].get_ylim()
	ax[4].set_ylim((0.,y4b))

	# Plot rtr
	ax[5].set_xscale('log')
	ax[5].plot(t[1], rtr[1]/pB[1], 'r-')
	if not any(np.isnan(rtr[2])):
		ax[5].plot(t[2], rtr[2]/pB[2], 'b-')
		ax[5].set_ylabel('Truncation Radius/pB)')
	(y5a, y5b) = ax[5].get_ylim()
	ax[5].set_ylim((0.,y5b))

	ax[1].yaxis.tick_right()
	ax[1].yaxis.set_ticks_position('both')
	ax[1].yaxis.set_label_position("right")
	ax[3].yaxis.tick_right()
	ax[3].yaxis.set_ticks_position('both')
	ax[3].yaxis.set_label_position("right")

	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	plt.savefig('TimeData/AvsT-'+d+'.png')
#	plt.savefig('TimeData/AvsT-'+d+'.eps',format='eps',dpi=1000)
	plt.clf()
	plt.close(f)

sumfile = open('TimeData/DiskSummary.txt','w')
sumfile.write(hdr+'\n')
for line in summary:
	sumfile.write(fmt.format(line)+'\n')
sumfile.close()
#ascii.write(summary,'TimeData/PeriSummary.txt')

### Array sof data I want stats on
#minp  = summary['minpC']
#pf    = summary['minpC']
#minpO = minp[summary['vers']==0]
#minp2 = minp[summary['vers']==1]
#minp3 = minp[summary['vers']==2]
#rtr_a = summary['rtr']
#rtr_2 = rtr_a[summary['vers']==1]
#rtr_3 = rtr_a[summary['vers']==2]
## rtr/pericenter
#rp_a  = rtr_a/minp


#stats = np.array( [ [np.median(l) for l in [minp, minpO, minp3]],
#					[   mad(l)    for l in [minp, minpO, minp3]],
#					[np.mean(l)   for l in [minp, minpO, minp3]],
#					[np.std(l)    for l in [minp, minpO, minp3]],
#					[np.min(l)    for l in [minp, minpO, minp3]],
#					[np.max(l)    for l in [minp, minpO, minp3]] ])
#print('     all minp, minp-O, minp-3')
#print(np.concatenate( ( np.array([
#	['med'],['mad'],['mn'],['std'],['min'],['max']], dtype='S6'), 
#	stats), axis=1) )
### Plot pericenters

#f2, ax2 = plt.subplots(2,2)

#ax2[0,0].set_xscale('log')
##ax2[l].set_yscale('log')
#ax2[0,0].scatter( 	[line.minpC for line in summary][1:],
#					[line.iMf   for line in summary][1:] )
#ax2[0,1].hist([line.minpC for line in summary][1:])
#ax2[0,1].set_xlabel('All min(pericenter)')
#ax2[1,0].hist(minpO)
#ax2[1,0].set_xlabel('Orig min(pericenter)')
#ax2[1,1].hist(minp3)
#ax2[1,1].set_xlabel('Trip min(pericenter)')

#plt.savefig('TimeData/Summaries.png')
##plt.savefig('TimeData/Summaries.eps',format='eps',dpi=1000)
#plt.clf()
#plt.close(f2)




