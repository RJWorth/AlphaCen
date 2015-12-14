#from cgs_constants import mSun,mEarth,mMoon,mMars,AU,day
#from cgs_constants import G as G_cgs
from mks_constants import mSun,mEarth,mMoon,mMars,AU,day
from mks_constants import deg2rad,rad2deg
from mks_constants import G as G_mks
import Disks as D
import Merc as M
import random as R
import numpy as np
from numpy import pi, sqrt, sin, cos, arccos, pi

### tolerance for floats to count as 'equal'
tol = 1.e-15

###############################################################################
class Obj(object):
	'''General object class, but actually use the subclasses, which can be
	either Asteroidal or Cartesian.'''

#-----------------------------------------------------------------------------#
	def __init__(self, name, mass, mCent, density, s):

		self.name    = name
		self.mass    = mass
		self.mCent   = mCent
		self.density = density
		self.s       = s

#-----------------------------------------------------------------------------#
	def peri(self):
		return(self.a()*(1-self.e()))
	def apo(self):
		return(self.a()*(1+self.e()))
#-----------------------------------------------------------------------------#
	def mr(self):
		'''Reduced mass: 1/mr = sum(1/m_i). Units unchanged.'''
		m = np.array([self.mass, self.mCent])
		mr   = 1./sum( [1./i for i in m] )
		return(mr)

#-----------------------------------------------------------------------------#
	def gm(self):
		'''Gravitational parameter: gm = G_mks * (mass+mCent). Takes mSun,
		returns	mks units.'''
		mTot = mSun*(self.mass+self.mCent)
		gm = G_mks*mTot
		return(gm)

#-----------------------------------------------------------------------------#
	def eps(self):
		"""Specific orbital energy of Body around Central object"""
		eps = -self.gm()/(2.*self.a()*AU)
		return eps

#-----------------------------------------------------------------------------#
	def P(self):
		"""Orbital period in days"""
		P = 2*pi * sqrt( (self.a()*AU)**3. / (G_mks*(self.mCent+self.mass)*mSun) )/day
		return P
#-----------------------------------------------------------------------------#
	def RH(self):
		"""Hill radius of Body"""
		RH = self.a() * (1-self.e()) * (self.mass/(3.*self.mCent))**(1./3.)
		return RH		

#-----------------------------------------------------------------------------#
	def RH2(self, other):
		"""Mutual Hill radius of two Bodies"""
  
		assert (self.mCent == other.mCent), 'Bodies do not have same central object'
		RH2 = ((self.mass+other.mass)/(3.*self.mCent))**(1./3.) * \
													((self.a()+other.a())/2.)
		return RH2		
#-----------------------------------------------------------------------------#
	def dr(self,other):
		'''Distance between two objects.'''

		dr = abs(self.a() - other.a())

		return(dr)

#-----------------------------------------------------------------------------#
	def IsClose(self,other,rh=10.):
		'''Check if two adjacent objects are within their mutual Hill radius
		of each other.'''


		dr = self.dr(other)
		RH = self.RH2(other)

		if (rh*RH >= dr):
			coll = True
		else:
			coll = False

		return coll

###############################################################################
class AsteroidalObj(Obj):
	'''Object in a mercury simulation, input in asteroidal coords.
	Units: AU (a), degrees (i, g, n, m), mSun(mass, mCen), g/cm^3 (density)'''

#-----------------------------------------------------------------------------#
	def __init__(self, 
			name='Earth', mass=mEarth/mSun, mCent = 1., density=5.51, 
			s=[0., 0., 0.],
			elements = [1.,0.,0., 0.,0.,0.]):
		super(AsteroidalObj,self).__init__(name, mass, mCent, density, s)
		self.elements = elements

#-----------------------------------------------------------------------------#
	def a(self):
		return(self.elements[0])
	def e(self):
		return(self.elements[1])
	def i(self):
		return(self.elements[2])
	def g(self):
		return(self.elements[3])
	def n(self):
		return(self.elements[4])
	def m(self):
		return(self.elements[5])
#-----------------------------------------------------------------------------#
	def CalcCartesian(self):
		'''Calculates Cartesian coords from current asteroidal coords.'''

	# Note: must convert to/from mks, and must put angles in radians
		a, e, i, g, n, m = self.a(), self.e(), self.i(), \
						   self.g(), self.n(), self.m()
		x, y, z, vx, vy, vz = M.Merc_El2X(
			[a*(AU), e, i*pi/180., g*pi/180., n*pi/180., m*pi/180.],
			[self.mCent*mSun, self.mass*mSun])
		x, y, z = x/(AU), y/(AU), z/(AU)
		vx,vy,vz = vx*day/(AU), vy*day/(AU), vz*day/(AU)
		pos = [x,  y,  z]
		vel = [vx, vy, vz]
		return(pos, vel)


###############################################################################
class CartesianObj(Obj):
	'''Object in a mercury simulation, input in Cartesian coords.
	Units: AU (pos), AU/day(vel), mSun(mass, mCen), g/cm^3 (density)'''

#-----------------------------------------------------------------------------#
	def __init__(self, 
			name='Earth', mass=mEarth/mSun, mCent = 1., density=5.51, 
			s=[0., 0., 0.],
			pos=[1.,0.,0.], vel=[0., 0.017204516775334737,0.]):
		super(CartesianObj,self).__init__(name, mass, mCent, density, s)
		### Vectors of position and velocity (3 letters = 3 numbers)
		self.pos = np.array(pos)
		self.vel = np.array(vel)

#-----------------------------------------------------------------------------#
### Magnitude of r and v (1 letter = 1 number)
	def r(self):
		return(M.mag( self.pos))
	def v(self):
		return(M.mag(self.vel))

#-----------------------------------------------------------------------------#
### Get individual coordinate values
	def x(self):
		return(self.pos[0])
	def y(self):
		return(self.pos[1])
	def z(self):
		return(self.pos[2])
	def vx(self):
		return(self.vel[0])
	def vy(self):
		return(self.vel[1])
	def vz(self):
		return(self.vel[2])
#-----------------------------------------------------------------------------#
	def dr(self,other):
		'''Distance between two objects.'''
		x1 = self.pos
		x2 = other.pos
		dr = M.mag(x1-x2)
		return(dr)

#-----------------------------------------------------------------------------#
	def CalcAsteroidal(self):
		'''Calculates asteroidal coords from current Cartesian coords.'''

		a = self.a()
		e = self.e()
		i = self.i()

		return( np.array([a,e,i]) )

#-----------------------------------------------------------------------------#
#       Asteroidal elements:
#-----------------------------------------------------------------------------#
	def a(self):
		'''Calculate semimajor axis of an object's orbit. Input mks, output AU'''
		a = -self.gm()/(2.*self.eps())/AU
		return(np.array(a))

#-----------------------------------------------------------------------------#
	def e(self):
		'''Calculate eccentricity of an object's orbit'''
		
		val = 1.+2.*self.eps()*M.mag(self.h())**2./self.gm()**2.
		if ((val < 0) & (-val < tol)):
			val = 0
		print(val)
		
		e = sqrt(val)
		return(np.array(e))

#-----------------------------------------------------------------------------#
	def i(self):
		'''Calculate inclination of an object's orbit. Takes mks, 
		returns degrees'''
		i = arccos(self.h()[2]/M.mag(self.h()))*180./pi
		return(np.array(i))

#-----------------------------------------------------------------------------#
#       Internal functions to calculate asteroidal elements:
#-----------------------------------------------------------------------------#
	def h(self):
		'''Calculate angular momentum. Takes AU/day, returns mks.'''
		h = np.cross(self.pos*AU,self.vel*AU/day)
		return(np.array(h))

#-----------------------------------------------------------------------------#
	def eps(self):
		'''Calculate specific orbital energy
		where r is the distance between two stars, v is the relative velocity.
		Final units should be mks (m^2/s^2).'''
		eps = (self.v()*AU/day)**2./2.-self.gm()/(self.r()*AU)
		return(eps)

###############################################################################
def mag(x):
	'''Takes vector, returns magnitude of that vector.'''

	xbar = sqrt(sum([xi**2. for xi in x]))
	return(xbar)

###############################################################################
def ReadInObjList(whichdir='In/',fname='big.in'): 
	'''Read in object list from big.in-type file'''

	infile=open(whichdir+fname,'r')
	AllLines=np.array(infile.readlines())
	infile.close()
	header = np.array([False for i in range(len(AllLines))])
	for ind,line in enumerate(AllLines):
		if line.startswith(')') | ('style' in line) | ('epoch' in line):
			header[ind] = True
	nh = sum(header)
	assert (sum(header==False)%4)==0
	nObj=int(len(AllLines[header==False])/4)
	style = AllLines[np.array([('style' in line) for line in AllLines])][0].split()[-1]

	objlist = []
	if style == 'Cartesian':
		for i in range(nObj):
			num = nh+i*4
			tag = np.array( AllLines[num].split())
			print(tag)
			name = tag[0]
			mass = float(tag[np.array(['m=' in word for word in tag])][0].strip('m='))
			dens = float(tag[np.array(['d=' in word for word in tag])][0].strip('d='))
			x = np.array(AllLines[num+1].split()).astype(np.float)
			v = np.array(AllLines[num+2].split()).astype(np.float)
			s = np.array(AllLines[num+3].split()).astype(np.float)
			objlist.append(M.CartesianObj(name=name, mass=mass, density=dens,
			pos=x,vel=v,s=s) )
	elif style == 'Asteroidal':
		for i in range(nObj):
			num = nh+i*4
			tag = np.array( AllLines[num].split())
			name = tag[0]
			mass = float(tag[np.array(['m=' in word for word in tag])][0].strip('m='))
			dens = float(tag[np.array(['d=' in word for word in tag])][0].strip('d='))
			a, e, i = np.array(AllLines[num+1].split()).astype(np.float)
			g, n, m = np.array(AllLines[num+2].split()).astype(np.float)
			elements = [a,e,i,g,n,m]
			s = np.array(AllLines[num+3].split()).astype(np.float)
			objlist.append(M.AsteroidalObj(name=name, mass=mass, density=dens,
			elements=elements,s=s) )

	return objlist

###############################################################################
def ReadObjData(name='Mars',whichdir='In/',fname='big.in'):
	'''Returns object containing the parameters for an object on a big.in or 
	small.in list'''

	objlist = M.ReadInObjList(whichdir=whichdir,fname=fname)
	indlist = []
	for ind, obj in enumerate(objlist):
		if obj.name==name:
			thisobj = obj
			indlist.append(ind)

	assert len(indlist)>0, 'No matching object!'
	assert len(indlist)<2, 'Multiple matching objects! indices: {}'.format(indlist) 

	return(thisobj)

###############################################################################
def MakeEjecShellObjList(whichdir='In/',cent='Mars',vmin=0.01, vmax=10., r=1.01): 
	'''Generate list of objects being ejected, in shell around obj cent,
	with v_inf ranging from vmin to vmax, at distance r in Hill radii'''

### Get coordinates of central object (assumes it is in big.in)
#	xcent, vcent, scent, mcent
	CenterOjb = M.ReadObjData(whichdir=whichdir,name=cent)

### Create positions of objects at distance r*rH from center of cent, 
### at random angles, with v_inf pulled randomly from [vmin,vmax]
	if type(CenterObj)=='CartesianObj':
		CenterObj.pos, CenterObj.vel
	elif type(CenterObj)=='AsteroidalObj':
		CenterObj.CalcCartesian()
	else:
		raise ValueError('Error: CenterObj type not recognized')

#	a, mass = disk.debris.ListParams()

	digits = str(len(str(len(a))))
	fmt = 'P{0:0'+digits+'}'

### List of planetesimals on circular co-planar orbits with random phases
	objlist = []
	for i in range(len(a)):
		objlist.append( M.Obj(name=fmt.format(i), mass=mass[i], density=3.,
		a=a[i], 
		g=R.uniform(0.,360.), n=R.uniform(0.,360.), m=R.uniform(0.,360.)) )

	return objlist

###############################################################################
def MakeDiskObjList(rtr = 5., sigC = 10., rh = 10., m = [mMoon/mSun,mMars/mSun], 
	f= [0.5, 0.5], alpha = 1.5, starA=True, iMax=180.):
	'''Generate list of objects based on the Disks semi-analytic model'''

	disk = D.Disk(r_out=rtr, alpha=alpha, sigC=sigC, rh=rh)
	disk.DebrisGenM(m,f)
	a, mass = disk.debris.ListParams()

#	digits = str(int(np.floor(np.log10(len(a))+1)))
	digits = str(len(str(len(a))))
	fmt = 'P{0:0'+digits+'}'

### List of planetesimals on circular co-planar orbits with random phases
	objlist = []
	for i in range(len(a)):
		objlist.append( M.Obj(name=fmt.format(i), mass=mass[i], density=3.,
		a=a[i], 
		g=R.uniform(0.,360.), n=R.uniform(0.,360.), m=R.uniform(0.,360.)) )

	### Add binary star (AlCenA)
	if starA==True:
		objlist.append( M.Obj(name='AlCenA', mass=1.105, density=1.5,
		a=23.7,               e=0.5179,             i=R.uniform(0.,iMax), 
		g=R.uniform(0.,360.), n=R.uniform(0.,360.), m=R.uniform(0.,360.)) )

	return objlist

###############################################################################
def WriteObjInFile(objlist='default', loc = 'Merc95/In/',infile='big', epoch=0.):
	'''Write a big.in or small.in file for mercury'''

### Currently non-variable variable
	style = 'Asteroidal'

### Make list of object parameters, unless provided
	if (objlist=='default'):
		objlist = M.MakeDiskObjList()

### Process big/small differences
	assert ((infile == 'big') | (infile=='small')), 'invalid infile: must be "big" or "small"'
	fname = loc+infile+'.in'
	if (infile == 'big'):
		vers = 0
	elif (infile == 'small'):
		vers = 1

	header = ')O+_06 '+['Big','Small'][vers]+ \
"-body initial data  (WARNING: Do not delete this line!!)\n) Lines beginning with `)' are ignored.\n)---------------------------------------------------------------------\n style (Cartesian, Asteroidal, Cometary) = {style}"\
+["\n epoch (in days) = {epoch}",""][vers]+ \
"\n)---------------------------------------------------------------------\n"

	objstr = '''  {0.name:16}  m={0.mass}  d={0.density}
    {0.a: .18e} {0.e: .18e} {0.i: .18e}
    {0.g: .18e} {0.n: .18e} {0.m: .18e}
    {0.s[0]: 19.18e} {0.s[1]: 19.18e} {0.s[2]: 19.18e}
'''

	with open(fname, 'w') as f:
		f.write(header.format(style=style,epoch=epoch))
		for i in range(len(objlist)):
			f.write(objstr.format(objlist[i]))

###############################################################################
def WriteParamInFile(loc = 'Merc95/In/', f = 'in', alg='hybrid', 
	ti=0., tf=365.25e3, tOut = 'default', dt=1., acc=1.e-12, 
	CEstop='no',CE='yes',CEfrag='no',tUnit='years',tRel='yes',prec='medium',
	rel='no',user='no',
	rEj=100, rStar=0.005, mStar=1.0, rCE=3.,dtDump='default',dtPer=100):
	'''Write a param.in file for mercury'''
	
	if (tOut == 'default'):
		tOut = tf/1.e3
	if (dtDump == 'default'):
		dtDump = int(round(tf/dt/100.))

	text = ''')O+_06 Integration parameters  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
) Important integration parameters:
)---------------------------------------------------------------------
 algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = {alg}
 start time (days)= {ti}
 stop time (days) = {tf}
 output interval (days) = {tOut}
 timestep (days) = {dt}
 accuracy parameter={acc}
)---------------------------------------------------------------------
) Integration options:
)---------------------------------------------------------------------
 stop integration after a close encounter = {CEstop}
 allow collisions to occur = {CE}
 include collisional fragmentation = {CEfrag}
 express time in days or years = {tUnit}
 express time relative to integration start time = {tRel}
 output precision = {prec}
 < not used at present >
 include relativity in integration= {rel}
 include user-defined force = {user}
)---------------------------------------------------------------------
) These parameters do not need to be adjusted often:
)---------------------------------------------------------------------
 ejection distance (AU)= {rEj}
 radius of central body (AU) = {rStar}
 central mass (solar) = {mStar}
 central J2 = 0
 central J4 = 0
 central J6 = 0
 < not used at present >
 < not used at present >
 Hybrid integrator changeover (Hill radii) = {rCE}
 number of timesteps between data dumps = {dtDump}
 number of timesteps between periodic effects = {dtPer}
'''.format(
	alg=alg, ti=ti, tf=tf, tOut=tOut, dt=dt, acc=acc,
	CEstop=CEstop, CE=CE, CEfrag=CEfrag, tUnit=tUnit, tRel=tRel,prec=prec,
	rel=rel, user=user,
	rEj=rEj, rStar=rStar,mStar=mStar, rCE=rCE, dtDump=dtDump, dtPer=dtPer)

	fname=loc+'param.'+f
	with open(fname, 'w') as f:
		f.write(text)


###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   units?
def Merc_El2X(el, mass):
	'''Convert orbital elements to cartesian for an ellipse (e < 1). 
Based on MCO_EL2X.FOR from J. Chambers' mercury6_2.for:
	gm = grav const * (central + secondary mass)
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of perihelion !!! )   in
	n = longitude of ascending node ) radians
	l = mean anomaly                )

	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(gm/a) )'''
# Big.in format:
# Asteroidal = Keplerian orbital elements, in an `asteroidal' format.
#              i.e.  a e I g n M, where
#               a = semi-major axis (in AU)
#               e = eccentricity
#               I = inclination (degrees)
#               g = argument of pericentre (degrees)
#               n = longitude of the ascending node (degrees)
#               m = mean anomaly (degrees)

### Extract needed parameters from input list
	a,e,i,g,n,m = [float(i) for i in el]
	gm = G_mks*sum(mass)

### Convert input degrees to radians
#	i, g, n, m = deg2rad*i, deg2rad*g, deg2rad*n, deg2rad*m

### Pericenter
	q = (1.-e)*a

### Get p (longitude of pericenter) from g (argument of pericenter)
	p = g + n

### Rotation factors
	si, ci = sin(i), cos(i)
	sg, cg = sin(g), cos(g)
	sn, cn = sin(n), cos(n)
	z1 = cg * cn
	z2 = cg * sn
	z3 = sg * cn
	z4 = sg * sn
	d11 =  z1 - z4*ci
	d12 =  z2 + z3*ci
	d13 =       sg*si
	d21 = -z3 - z2*ci
	d22 = -z4 + z1*ci
	d23 =       cg*si

	assert (e < 1.0)	# haven't finished the other cases yet
### Calculate ellipse
	if (e < 1.0):
	### Ellipse
		romes = sqrt(1.0 - e*e)
		temp = M.Merc_KeplerEllipse(e,m)
		se, ce = sin(temp), cos(temp)
		z1 = a * (ce - e)
		z2 = a * romes * se
		temp = sqrt(gm/a) / (1.0 - e*ce)
		z3 = -se * temp
		z4 = romes * ce * temp
#	elif (e == 1.0) then
	### Parabola
#		ce = orbel_zget(l)
#		z1 = q * (1.0 - ce*ce)
#		z2 = 2.0 * q * ce
#		z4 = sqrt(2.0*gm/q) / (1.0 + ce*ce)
#	 	z3 = -ce * z4
#	elif (e > 1.):
	### Hyperbola
#		romes = sqrt(e*e - 1.0)
#		temp = orbel_fhybrid(e,l)
#		se, ce = sin(temp), cos(temp)
#		z1 = a * (ce - e)
#		z2 = -a * romes * se
#		temp = sqrt(gm/abs(a)) / (e*ce - 1.0)
#		z3 = -se * temp
#		z4 = romes * ce * temp
	
###
	x = d11 * z1  +  d21 * z2
	y = d12 * z1  +  d22 * z2
	z = d13 * z1  +  d23 * z2
	u = d11 * z3  +  d21 * z4
	v = d12 * z3  +  d22 * z4
	w = d13 * z3  +  d23 * z4

	return(x,y,z,u,v,w)

###############################################################################
def Merc_KeplerEllipse(e,oldl):
	'''Solves Kepler's equation for eccentricities less than one.
 Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.

  e = eccentricity
  l = mean anomaly      (radians)
  u = eccentric anomaly (   "   )'''
#------------------------------------------------------------------------------

	from math import pi, exp, log, sin, cos, sqrt
	twopi = 2.*pi
	piby2 = pi/2.

### Reduce mean anomaly to lie in the range 0 < l < pi
	if (oldl >= 0):
		l = oldl % twopi
	else:
		l = oldl % twopi + twopi
	sign = 1.0
	if (l > pi):
		l = twopi - l
		sign = -1.0

	ome = 1.0 - e
#-----------------------------
	if ((l >= .45) | (e < .55)):
### Regions A,B or C in Nijenhuis

### Rough starting value for eccentric anomaly
		if (l < ome):
			u1 = ome
		elif (l > (pi-1.0-e)):
			u1 = (l+e*pi)/(1.0+e)
		else:
			u1 = l + e

### Improved value using Halley's method
		flag = u1 > piby2
		if (flag):
			x = pi - u1
		else:
			x = u1
		x2 = x*x
		sn = x*(1.0 + x2*(-.16605 + x2*.00761) )
		dsn = 1.0 + x2*(-.49815 + x2*.03805)
		if (flag):
			dsn = -dsn
		f2 = e*sn
		f0 = u1 - f2 - l
		f1 = 1.0 - e*dsn
		u2 = u1 - f0/(f1 - .5*f0*f2/f1)

#---------------------
	else:
### Region D in Nijenhuis
### Rough starting value for eccentric anomaly
		z1 = 4.0*e + .50
		p = ome / z1
		q = .50 * l / z1
		p2 = p*p
		z2 = exp( log( sqrt( p2*p + q*q ) + q )/1.5 )
		u1 = 2.0*q / ( z2 + p + p2/z2 )

### Improved value using Newton's method
		z2 = u1*u1
		z3 = z2*z2
		u2 = u1 - .075*u1*z3 / (ome + z1*z2 + .375*z3)
		u2 = l + e*u2*( 3.0 - 4.0*u2*u2 )

### Accurate value using 3rd-order version of Newton's method
### N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!

### First get accurate values for u2 - sin(u2) and 1 - cos(u2)
	bigg = (u2 > piby2)
	if (bigg):
		z3 = pi - u2
	else:
		z3 = u2

	big = (z3 > (.5*piby2))
	if (big):
		x = piby2 - z3
	else:
		x = z3

	x2 = x*x
	ss = 1.0
	cc = 1.0

	ss = x*x2/6.*(1. - x2/20. *(1. - x2/42. *(1. - x2/72. *(1. -
         x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
	cc =   x2/2.*(1. - x2/12. *(1. - x2/30. *(1. - x2/56. *(1. -
         x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
         x2/306.))))))))

	if (big):
		z1 = cc + z3 - 1.0
		z2 = ss + z3 + 1.0 - piby2
	else:
		z1 = ss
		z2 = cc

	if (bigg):
		z1 = 2.0*u2 + z1 - pi
		z2 = 2.0 - z2

	f0 = l - u2*ome - e*z1
	f1 = ome + e*z2
	f2 = .5*e*(u2-z1)
	f3 = e/6.0*(1.0-z2)
	z1 = f0/f1
	z2 = f0/(f2*z1+f1)
	mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )

	return(mco_kep)

###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   for orbit with no inclination
def El2X(el, mass):
	'''Convert orbital elements to cartesian for an ellipse (e < 1) with zero inclination
	gm = grav const * sum of masses
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of pericenter     )   in
	la = longitude of ascending node) radians
	f = true anomaly                )
	E = eccentric anomaly           )
	m = mean anomaly                )


	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
	EDIT: USE MKS UNITS TO MATCH G'''


### Extract needed parameters from input list
	a,e,i,g,la,f = [float(i) for i in el]
	gm = G_mks*sum(mass)

### Convert input degrees to radians
#	i, g, la = deg2rad*i, deg2rad*g, deg2rad*la
	f = deg2rad*f

### Calculate radius
	r = a*(1-e**2)/(1+e*cos(f))

### Position coords
	x = r*cos(f)
	y = r*sin(f)
	z = 0.

### Period
	T = ( (4*pi**2/gm) * a**3 )**0.5

### Mean motion
	n = 2*pi/T

### Velocity coords
	u = -    sin(f)  * n*a/(1-e**2)**0.5
	v =   (e+cos(f)) * n*a/(1-e**2)**0.5
	w = 0.

	return(x,y,z,u,v,w)
	
###########################################################################
def WriteRuntime(t1,t2):
	'''Return the sim runtime formatted readably'''

	secs = t2-t1
	mins = secs/60.
	hrs  = mins/60.
	days = hrs/24.

	if (mins <= 1.):
		return( '--------- runtime = {0} {1} ---------'.format(secs,'sec') )
	elif (hrs <= 1.):
		return( '--------- runtime = {0:.2f} {1} ---------'.format(mins,'min') )
	elif (days <= 1.):
		return( '--------- runtime = {0:.2f} {1} ---------'.format(hrs, 'hrs') )
	else:
		return( '--------- runtime = {0:.2f} {1} ---------'.format(days, 'days') )

###########################################################################
