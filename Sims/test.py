

import AlphaCenModule
import numpy as np
from numpy import log10, sqrt
AU		= 1.49597871e11	# m
mSun		= 1.9891e30		# kg
mA, mB, mC	= 1.105*mSun, 0.934*mSun, 0.123*mSun
m			= (mA, mB, mC)
ind1, ind2 = -5, -1
nsteps=ind2-ind1
xvA = np.array([[0.,0.,0., 0.,0.,0.] for i in range(nsteps)])
xvB = AlphaCenModule.ReadAei('Prx3/6-3', 'AlCenB', ind1, ind2)
xvC = AlphaCenModule.ReadAei('Prx3/6-3', 'PrxCen', ind1, ind2)

xvCM, xvAcm, xvBcm, xvCcm = AlphaCenModule.CM(xvA, xvB, xvC, m)

rA, EA, rB, EB, rC, EC, rAB = AlphaCenModule.GetRE(xvAcm, xvBcm, xvCcm, m)

reload(AlphaCenModule)
mu = AlphaCenModule.Mu(m)











