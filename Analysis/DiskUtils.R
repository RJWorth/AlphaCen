### Get energies
G = 6.67384e-11		# m^3/(kg*s^2)
mSun = 1.9891e30	# kg
m1 = 0.123*mSun
m2 = 1000.
mu = G*(m1+m2)

day = 3600*24.		# s
AU  = 149597870700 	# m

###############################################################################
sum3d = function(matrix3d)	{
	flattened = matrix3d[1,,]

	dims = dim(matrix3d)
	for (i in 2:dims[1])	{
		flattened = rbind(flattened,matrix3d[i,,])
	}
	
	summ = summary(flattened)
	return(summ)
	}


r = function(i,j)	{
	r = sqrt(disk[i,j, 'x']^2+disk[i,j, 'y']^2)*AU
	return(r)
	}

v = function(i,j)	{
	v = sqrt(disk[i,j,'vx']^2+disk[i,j,'vy']^2)*AU/day
	return(v)
	}

keps = function(i,j)	{
	keps = (v(i,j)^2)/2.
	return(keps)
	}
veps = function(i,j)	{
	veps = - mu/r(i,j)
	return(veps)
	}
eps = function(i,j)	{
	eps = keps(i,j)+veps(i,j)
	return(eps)
	}

enorm = function(i,j) {
	enorm = eps(i,j)/mean(eps(i,1:1001))
	return(enorm)
	}

vorb_expect = function(i,j) {
	vorb = sqrt( mu/r(i,j))
	return(vorb)
	}

a = function(i,j) {
	a = -mu/(2.*eps(i,j))
	return(a)
	}
avar = function(i,j) {
	a = a(i,j)
	avar = (max(a)-min(a))/mean(a)
	return(avar)
	}


