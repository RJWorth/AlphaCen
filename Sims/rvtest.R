
dir='TestC'
aeidir=paste(dir,'/Out/AeiOutFiles/',sep='')

stars=list()
starnames=c('AlCenB')
for (j in 1:length(starnames)) stars[[j]]=read.table(
	paste(aeidir,starnames[j],'.aei',sep=''), header=F,skip=4,
	col.names=c('Time','a','e','i','mass','dens', 'x','y','z','vx','vy','vz')
	)[998:1001,c(1:3,7:12)]
t=stars[[1]]$Time

attach(stars[[1]])

### Constants
G	= 6.67384e-11			# m^3/kg/s^2
AU	= 1.49597871e11			# m
day	= 24*3600.				# seconds
Msun= 1.9891e30				# kg

### Calculated values
m=Msun*c(1.105, .934)
mu=G*sum(m)

# r, v in AU units
r0=sqrt( x^2 + y^2 + z^2)	# AU
v0=sqrt(vx^2 +vy^2 +vz^2)	# AU/day

# calculate a from v/r:
a0exp = 1/( 2/r0 - v0^2/(mu*day^2/AU^3) )
v0exp=sqrt((mu*day^2/AU^3)*(2/r0 - 1/mean(a0exp)))

# r, v in mks
r1=r0*AU					# m
v1=v0*AU/day				# m/s

# calculate a from v/r:
aexp = 1/( 2/r1 - v1^2/mu )
# expected v(r)
vexp=sqrt(mu*(2/r1 - 1/mean(aexp)))

pdf(paste(dir,'/Out/AeiOutFiles/test1.pdf',sep=''), height=4, width=4)
plot(r0,v0,type='p',pch=20)
points(r0,v0exp,type='p',pch='.',col='blue')
dev.off()


pdf(paste(dir,'/Out/AeiOutFiles/test2.pdf',sep=''), height=4, width=4)
plot(r1,v1,type='p',pch=20)
points(r1,vexp,type='p',pch='.',col='blue')
dev.off()


detach(stars[[1]])



