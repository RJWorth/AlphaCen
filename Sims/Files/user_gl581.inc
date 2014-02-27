c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      USER.INC    (ErikSoft   4 March 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Avi Mandell
c
c Parameters that govern migration:
c
c   NMIG = # of migrating planets, in order in BIG file
c   TMAX = Time that migration ends, in days
c   MIGFACT = Ratio of  migration time to TMAX for each planet; 
c       if two planets have MIGFACT in the same ratio as their distances, 
c       they will migrate with the same inward velocity
c   TMIG0 = Init. mig time, as a fract. of total time
c   SWITCH = Time to switch to expon. mig fade-out, as a fract. of total time
c   GAMMA = Var. mig exponent
c   STOPDIST = Dist. to stop mig at, if reached
c
c   Gas Drag Formulas From Thommes03:
c   GASTIME = Time at which gas drag ceases, in days
c   RHOCONST = Gas density const, in g/cm^3
c   RHOEXP = Gas profile expon.
c   ZCONST = Vert. gas density scale height const., in AU
c   ZEXP = Vert. gas profile expon.
c   CDRAG = Drag constant
c   MEFF = Effective part. mass, in g
c   GCONST = Gravitational const. in AU,Sol.masses,days
c
      integer nmig
      real*8 tmig0,tmax,migfact(1),switch,gamma,stopdist
c
      real*8 gastime,rhoconst,rhoexp,zconst,zexp
      real*8 cdrag,meff,gconst
c
      parameter (nmig = 0)                          
      parameter (tmax = 36525000.00)
      data migfact /1.00d0/		       
      parameter (tmig0 = 0.00d0)		!Stan.:0.67         
      parameter (switch = 0.5d0)    		!Stan.:0.5
      parameter (gamma = 0.6d0)			!Stan.:0.6
      parameter (stopdist = 0.045d0)                 
						    
      parameter (gastime = 0)
      parameter (rhoconst = 1.4d-9)            	!Stan.:1.4d-9
      parameter (rhoexp = (11.0/4.0))              
      parameter (zconst = 0.0472d0)                 
      parameter (zexp = (5.0/4.0))                  
      parameter (cdrag = 1.0d0)                     
      parameter (meff = 1.26d16)		!1e-23 M_e
      parameter (gconst = 2.9601e-4)
