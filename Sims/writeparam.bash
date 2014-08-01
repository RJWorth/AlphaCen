#!/bin/sh
###############################################################################

if [ $2 = $5 ]; then
	paramfile=$1/In/param.in
else
	paramfile=$1/Out/param.dmp
fi

echo '	writeparam.bash,          stop time = 365.25e'$2

dumpsteps=$(echo "365.25*10^$2/$4/5" | bc)
#printf "writeparam dumpsteps =  %20s \n" $dumpsteps >> testing.txt

output=$(echo "$2-3" | bc)

echo ')O+_06 Integration parameters  (WARNING: Do not delete this line!!)' >  $paramfile
echo ') Lines beginning with ")" are ignored.                   ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ') Important integration parameters:                       ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ' algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = hyb        ' >> $paramfile
echo ' start time (days)= 0.0                                   ' >> $paramfile
echo ' stop time (days) = 365.25e'$2'                           ' >> $paramfile
#echo ' output interval (days) = 365.25e'$3'                     ' >> $paramfile
echo ' output interval (days) = 365.25e'$output'                ' >> $paramfile
echo ' timestep (days) = '$4'                                   ' >> $paramfile
echo ' accuracy parameter=1.d-12                                ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ') Integration options:                                    ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ' stop integration after a close encounter = no            ' >> $paramfile
echo ' allow collisions to occur = yes                          ' >> $paramfile
echo ' include collisional fragmentation = no                   ' >> $paramfile
echo ' express time in days or years = years                    ' >> $paramfile
echo ' express time relative to integration start time = yes    ' >> $paramfile
echo ' output precision = med                                   ' >> $paramfile
echo ' < not used at present >                                  ' >> $paramfile
echo ' include relativity in integration= no                    ' >> $paramfile
echo ' include user-defined force = '$6'                        ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ') These parameters do not need to be adjusted often:      ' >> $paramfile
echo ')---------------------------------------------------------' >> $paramfile
echo ' ejection distance (AU)= 100000                           ' >> $paramfile
echo ' radius of central body (AU) = 0.0051                     ' >> $paramfile
echo ' central mass (solar) = '$7'                              ' >> $paramfile
echo ' central J2 = 0                                           ' >> $paramfile
echo ' central J4 = 0                                           ' >> $paramfile
echo ' central J6 = 0                                           ' >> $paramfile
echo ' < not used at present >                                  ' >> $paramfile
echo ' < not used at present >                                  ' >> $paramfile
echo ' Hybrid integrator changeover (Hill radii) = 3.           ' >> $paramfile
echo ' number of timesteps between data dumps = '$dumpsteps'    ' >> $paramfile
echo ' number of timesteps between periodic effects = 10000     ' >> $paramfile
