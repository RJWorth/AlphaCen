


### Define values for constants
AU			= 1.49597871e13	# cm
day			= 24*3600.		# s
G			= 6.67384e-8	# cm^3/g/s^2 
mSun		= 1.9891e33		# g
mEarth      = 5.97219e27    # g
mMoon       = 7.34767309e25 # g
mMars       = 6.41693e26    # g

### Masses for AlphaCen system
mA, mB, mC	= 1.105*mSun, 0.934*mSun, 0.123*mSun
m 			= [mA, mB, mC]
M 			= mA+mB+mC
mu 			= G*M


