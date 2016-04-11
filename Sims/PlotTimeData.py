
import numpy as np
import AlphaCenModule as AC
import Merc
import matplotlib.pyplot as plt
from astropy.io import ascii
from statsmodels.robust.scale import mad as madfn
import pickle
import pandas
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

### Laptop data locations
Dirs1 = ['d'+'{0:02}'.format(i) for i in range(1,35)]
Dirs2 = ['a'+   '{0}'.format(i) for i in range(1,13)]
Dirs3 = ['b'+'{0:02}'.format(i) for i in range(1,27)]
Dirs4 = ['c'+'{0:02}'.format(i) for i in range(1,9)]

Dirs = Dirs2+Dirs3+Dirs4+Dirs1
#Dirs = Dirs[0:2]
Dirs = ['d05','d10']
ndir = len(Dirs)

### Code for Original, DiskB-2, and DiskB-3 versions of each sim
vers = ['O','A2','A3','B2','B3']

### Object class to easily fill one line of data into summary format string
#-----------------------------------------------------------------------------#
#class Line(object):
#	'''Line of data on one sim/vers, for summary.'''
#	def __init__(self, dr,vrs,tf,
#				aBf,eBf, aCf,eCf,
#				minpB,mintB,minpC,mintC,
#				iMf,rtrf):
#		self.dr     = dr
#		self.vrs   = vrs
#		self.tf     = tf
#		self.aBf    = aBf
#		self.eBf    = eBf
#		self.aCf    = aCf
#		self.eCf    = eCf
#		self.minpB  = minpB
#		self.mintB  = mintB
#		self.minpC  = minpC
#		self.mintC  = mintC
#		self.iMf    = iMf
#		self.rtrf   = rtrf
#-----------------------------------------------------------------------------#
# names of columns in summary
hdrcols =  ['dr','vrs',  'tf',  'aBf', 'eBf', 'pBf',  'apBf',
								'aCf', 'eCf', 'pCf',  'apCf',
					'minpB','mintB','minpC','mintC', 'iMf','rtrf','rpf','rpm']
# dtype for each column in summary as structured array
#sumformat= [' >S3',' >S3','9.2e',' 6.2f',' 8.4f',' 6.2f',' 6.2f',
#				'8.1f','8.4f','8.1f','8.1f',
#					  '7.2f', '9.2e', '7.2f', '9.2e','6.1f','6.2f','7.4f']
sumdtype  = {'names':hdrcols, 
'formats':[ 'S3','S3', 'f', 'f','f','f','f',
				'f','f','f','f', 
			'f','f', 'f', 'f','f','f','f','f']}
# column widths (can this be better extracted from sumdtype?
widths  =  [   3,    3,     9,      6,      8,      6,      6,    
									8,     8,       8,      8,      
					7,      9,     7,   9, 6, 6, 7,7]
assert len(hdrcols) == len(widths) == len(sumdtype['formats']), \
	'lengths of format strings not equal! {0} {1} {2}'.format(
	   len(hdrcols) == len(widths) == len(sumdtype['formats']))
### Format string for columns wanted in summary
fmt=' '.join( ['{0.[''+hdrcols[i]+'']:'+sumdtype['formats'][i]+'}' for i in range(len(hdrcols)) ] )
# make header string of column names to print at top of file
hdrstr = ' '.join(['{0['+str(i)+']: >'+str(w)+'}' for i,w in enumerate(widths)])
hdr = hdrstr.format(hdrcols)
summary=np.zeros( (ndir*3), dtype = sumdtype)

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

	# calculate minimum pericenter, last time it happens, and final i
	for ind,l in enumerate([0,1,2]):
		# get basic stats of this sim
		dr    = d
		vrs   = vers[l]
		tf    = t[l][-1]
		aBf   = aB[l][-1]
		eBf   = eB[l][-1]
		pBf   = aBf*(1-eBf)
		apBf  = aBf*(1+eBf)
		aCf   = aC[l][-1]
		eCf   = eC[l][-1]
		pCf   = aCf*(1-eCf)
		apCf  = aCf*(1+eCf)
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
		rpf   = rtrf/pBf
		rpm   = rtrf/minpB
		# add to summary
#		line = np.array([[dr,vrs,tf, aBf,eBf,pBf,apBf,
#									 aCf,eCf,pCf,apCf,
#						minpB,max(mintB),minpC,max(mintC),iMf,rtrf]])
		summary[3*i+l] = dr,vrs,tf, aBf,eBf,pBf,apBf, aCf,eCf,pCf,apCf,\
						minpB,max(mintB),minpC,max(mintC),iMf,rtrf,rpf,rpm

	### Plot in eps, pdf, png
#	for itr in [1,2,3]:
	for itr in [2]:
		##### Make plots
		f, ax = plt.subplots(3, sharex=True, figsize=(5,5))
		### Plot star orbits for Original and DiskB-3
		axtwin = ['','','']
		for ind,l in enumerate([0]):
			print(d)
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
			ax[ind].set_ylabel('Distance (AU)')
		### Plot i over time
		ax[1].set_xscale('log')
		ax[1].plot(t[0], iB[0], 'r-')
		if not any(np.isnan(aC[0])):
			ax[1].plot(t[0], iM[0], 'k-',
						 t[0], iC[0], 'b-')
		ax[1].plot(t[2], iB[2], 'r--')
		if not any(np.isnan(aC[2])):
			ax[1].plot(t[2], iM[2], 'k--',
						 t[2], iC[2], 'b--')
		ax[1].set_xlabel('t (yr)')
		ax[1].set_ylabel('i (deg)')
		# Plot e
		ax[2].set_xscale('log')
		ax[2].plot(t[0], eB[0], 'r-')
		if not any(np.isnan(aC[0])):
			ax[2].plot(t[0], eC[0], 'b-')
		ax[2].plot(t[2], eB[2], 'r--')
		if not any(np.isnan(aC[2])):
			ax[2].plot(t[2], eC[2], 'b--')
		ax[2].set_ylabel('Eccentricity')

	#	# Plot rtr
	#	ax[4].set_xscale('log')
	#	ax[4].plot(t[1], rtr[1], 'r-')
	#	if not any(np.isnan(rtr[2])):
	#		ax[4].plot(t[2], rtr[2], 'b-')
	#		ax[4].set_ylabel('Truncation Radius (AU)')
	#	(y4a, y4b) = ax[4].get_ylim()
	#	ax[4].set_ylim((0.,y4b))

	#	# Plot rtr
	#	ax[5].set_xscale('log')
	#	ax[5].plot(t[1], rtr[1]/pB[1], 'r-')
	#	if not any(np.isnan(rtr[2])):
	#		ax[5].plot(t[2], rtr[2]/pB[2], 'b-')
	#		ax[5].set_ylabel('Truncation Radius/pB)')
	#	(y5a, y5b) = ax[5].get_ylim()
	#	ax[5].set_ylim((0.,y5b))

		ax[1].yaxis.tick_right()
		ax[1].yaxis.set_ticks_position('both')
		ax[1].yaxis.set_label_position("right")
	#	ax[3].yaxis.tick_right()
	#	ax[3].yaxis.set_ticks_position('both')
	#	ax[3].yaxis.set_label_position("right")

		ax[-1].set_xlabel('Time (yrs)')
	#	plt.gcf().subplots_adjust(bottom=0.15)
		f.subplots_adjust(hspace=0)
		plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
		if itr == 1:
			plt.savefig('../Paper/Inserts/AvsT-'+d+'.eps',format='eps',dpi=1000)
		elif itr == 2:
			plt.savefig('../Paper/Inserts/AvsT-'+d+'.pdf')
		elif itr == 3:
			plt.savefig('../Paper/Inserts/AvsT-'+d+'.png', dpi=1000)
		else:
			assert 0==1, 'invalid iteration, no AvsT plot saved!'
		plt.clf()
		plt.close(f)

### Convert summary to record array
summary = np.rec.array(summary)
from prettytable import PrettyTable
prettysum = PrettyTable(summary.dtype.names)
for row in summary:
	prettysum.add_row(row)
prettysum.align='r'
### Save summary in human-readable text file
if ndir <= 80:
	sumfile = open('TimeData/DiskSummary.txt','w')
	#sumfile.write(hdr+'\n')
	#for line in summary:
	#	sumfile.write(fmt.format(line)+'\n')
	sumfile.write(str(prettysum))
	sumfile.close()
	### Pick summary to read back into python
	pickle.dump( summary, open('TimeData/DiskSummary.pkl','wb') )

#ascii.write(summary,'TimeData/PeriSummary.txt')

### Arrays of data I want stats on
# masks
prx = np.array(summary['apCf']>10000.)
v0  = np.array(summary['vrs']=='O')
v2  = np.array(summary['vrs']=='2')
v3  = np.array(summary['vrs']=='3')
sml = np.array(summary['rpf']<0.2)

with open('TimeData/Stats.txt','w') as StatsFile:
	StatsFile.write('Truncation Radius Statistics for various slices\n')

#-----------------------------------------------------------------------------#
def PrintStats(x, f=''):
	if (len(x) ==0):
		assert 0==1, 'nothing passed to fn'
	x = np.array(x)
	nonanx = np.array([~np.isnan(x)])
	l    = len(x)
	nnan = sum(np.isnan(x))
	med  = np.median(x[~np.isnan(x)])
	mad  = np.median(abs( x[~np.isnan(x)] - med ))/0.6745
	mean = np.nanmean(x)
	std  = np.nanstd(x)
	xmin = np.nanmin(x)
	xmax = np.nanmax(x)
	if f=='':
		print('{l:3} objs, {n:3} nans ({xmin:5.3f}--{xmax:5.3f})\n'.format(
											l=l,n=nnan,xmin=xmin,xmax=xmax))
		print('med = {med:5.3f} +- {mad:5.3f}\n'.format(med=med,mad=mad))
		print('mn  = {mean:5.3f} +- {std:5.3f}\n'.format(mean=mean, std=std))
	else: 
		with open('TimeData/Stats.txt','a') as StatsFile:
			f.write('{l:3} objs, {n:3} nans ({xmin:5.3f}--{xmax:5.3f})\n'.format(
											l=l,n=nnan,xmin=xmin,xmax=xmax))
			f.write('med = {med:5.3f} +- {mad:5.3f}\n'.format(med=med,mad=mad))
			f.write('mn  = {mean:5.3f} +- {std:5.3f}\n'.format(mean=mean, std=std))

#-----------------------------------------------------------------------------#
with open('TimeData/Stats.txt','a') as StatsFile:
	StatsFile.write('-------------------rtrf-All-------------------------\n')
	PrintStats(summary[:]['rtrf'],StatsFile)
	StatsFile.write('-------------------rtrf-v2--------------------------\n')
	PrintStats(summary[v2]['rtrf'],StatsFile)
	StatsFile.write('-------------------rtrf-v3--------------------------\n')
	PrintStats(summary[v3]['rtrf'],StatsFile)
	StatsFile.write('-------------------rtrf-v3 & prx--------------------\n')
	PrintStats(summary[prx & v3]['rtrf'],StatsFile)

	StatsFile.write('\n-------------------rpf---All-------------------------\n')
	PrintStats(summary[:]['rpf'],StatsFile)
	StatsFile.write('-------------------rpf---v2--------------------------\n')
	PrintStats(summary[v2]['rpf'],StatsFile)
	StatsFile.write('-------------------rpf---v3--------------------------\n')
	PrintStats(summary[v3]['rpf'],StatsFile)
	StatsFile.write('-------------------rpf---v3 & prx--------------------\n')
	PrintStats(summary[prx & v3]['rpf'],StatsFile)

	StatsFile.write('\n-------------------rpm---All-------------------------\n')
	PrintStats(summary[:]['rpm'],StatsFile)
	StatsFile.write('-------------------rpm---v2--------------------------\n')
	PrintStats(summary[v2]['rpm'],StatsFile)
	StatsFile.write('-------------------rpm---v3--------------------------\n')
	PrintStats(summary[v3]['rpm'],StatsFile)
	StatsFile.write('-------------------rpm---v3 & prx--------------------\n')
	PrintStats(summary[prx & v3]['rpm'],StatsFile)

	StatsFile.write('\n-------------------iM---All-------------------------\n')
	PrintStats(summary[:]['iMf'],StatsFile)
	StatsFile.write('-------------------iM---Sml-------------------------\n')
	PrintStats(summary[sml]['iMf'],StatsFile)
	StatsFile.write('-------------------iM---Nrm-------------------------\n')
	PrintStats(summary[~sml]['iMf'],StatsFile)


#fmt=' '.join( ['{0.[''+hdrcols[i]+'']:'+sumdtype['formats'][i]+'}' for i in range(len(hdrcols)) ] )
#summary=np.zeros( (ndir*3), dtype = sumdtype)



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

### Plot summary data
nsamp = 50
rV2Prx    = summary[v2 &  prx]['rtrf']
rV2NonPrx = summary[v2 & ~prx]['rtrf']
rV3Prx    = summary[v3 &  prx]['rtrf']
rV3NonPrx = summary[v3 & ~prx]['rtrf']
rbins = np.linspace(np.nanmin(summary['rtrf']), np.nanmax(summary['rtrf']), nsamp+1)

rpV2Prx    = summary[v2 &  prx]['rpf']
rpV2NonPrx = summary[v2 & ~prx]['rpf']
rpV3Prx    = summary[v3 &  prx]['rpf']
rpV3NonPrx = summary[v3 & ~prx]['rpf']

rpbins = np.linspace(np.nanmin(summary['rpf']), np.nanmax(summary['rpf']), nsamp+1)


#f2, ax2 = plt.subplots(3,2)

## histograms of rtrf (AU) rp (r/pB)
#if len(rV2NonPrx) != 0:
#	ax2[0,0].hist( rV2NonPrx, rbins)
#	ax2[0,0].legend()
#	ax2[0,0].set_title('rTr: B-2')
#if len(rV3Prx) != 0:
#	ax2[1,0].hist( rV3Prx, rbins)
#	ax2[1,0].legend()
#	ax2[1,0].set_title('rTr: B-3, Prx')
#if len(rV3NonPrx) != 0:
#	ax2[2,0].hist( rV3NonPrx, rbins)
#	ax2[2,0].legend()
#	ax2[2,0].set_title('rTr: B-3, Non')

#(x1a, x2a) = ax2[0,0].get_xlim()
#(x1b, x2b) = ax2[1,0].get_xlim()
#(x1c, x2c) = ax2[2,0].get_xlim()
#ax2[0,0].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )
#ax2[1,0].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )
#ax2[2,0].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )

#if len(rpV2NonPrx) != 0:
#	ax2[0,1].hist( rpV2NonPrx, rpbins)
#	ax2[0,1].legend()
#	ax2[0,1].set_title('rP: B-2')
#if len(rpV3Prx) != 0:
#	ax2[1,1].hist( rpV3Prx, rpbins)
#	ax2[1,1].legend()
#	ax2[1,1].set_title('rP: B-3, Prx')
#if len(rpV3NonPrx) != 0:
#	ax2[2,1].hist( rpV3NonPrx, rpbins)
#	ax2[2,1].legend()
#	ax2[2,1].set_title('rP: B-3, Non')

#(x1a, x2a) = ax2[0,1].get_xlim()
#(x1b, x2b) = ax2[1,1].get_xlim()
#(x1c, x2c) = ax2[2,1].get_xlim()
#ax2[0,1].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )
#ax2[1,1].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )
#ax2[2,1].set_xlim( (min(x1a,x1b,x1c), max(x2a,x2b,x1c)) )

#plt.savefig('TimeData/Summaries.png')
##plt.savefig('TimeData/Summaries.eps',format='eps',dpi=1000)
#plt.clf()
#plt.close(f2)




