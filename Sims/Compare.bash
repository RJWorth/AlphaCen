#!/bin/bash
###############################################################################
# compare the orbital parameters in the original and disk triple systems

#module load python/2.7.3
#module load gcc/4.7.1
#module load R

h=$(pwd)
echo $h

hom='Proxlike/'
pre='Prx'
Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
##Dir=(01 03 05 06 07 09 10 11 12 13 14 15 16 17 19 20 21 22 23 24 26 33)
##Dir=(02 04 08 18 25 27 28 29 30 31 32 34)

#hom='Proxlike/071714/'
#pre='Prx'
#Dir=(1 2 3 4 5 6 7 8 9 10 11 12)

#hom='Proxlike/072314/'
#pre='Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)

#hom='Proxlike/081414/'
#pre='Prx'
#Dir=(01 02 03 04 05 06 07 08)

subdirs=(/DiskA-2 /DiskA-3 /DiskB-2 /DiskB-3)

for i in ${Dir[*]}; do
	d=$hom$pre$i
	echo '--------------------------------'$d'--------------------------------'

	for j in ${subdirs[*]}; do
	mkdir $d$j
	mkdir $d$j/In
	mkdir $d$j/Out
	scp -rp rjw274@hammer10.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d$j'/In/big.in'       $d$j'/In'
	scp -rp rjw274@hammer10.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d$j'/In/small.in'     $d$j'/In'
	scp -rp rjw274@hammer10.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d$j'/Out/info.out'    $d$j'/Out'
	scp -rp rjw274@hammer10.rcc.psu.edu:'~/work/AlphaCen/Sims/'$d$j'/Out/AeiOutFiles' $d$j'/Out'
	done

#	cd $d/Original/Out
#	scp -rp rjw274@shapiro.astro.psu.edu:'~/AlphaCen/Sims/'$d'/Original/Out/AeiOutFiles' .
#	scp -p  rjw274@shapiro.astro.psu.edu:'~/AlphaCen/Sims/'$d'/Original/Out/info.out' .
#	./elem
#	mv *.aei AeiOutFiles
#	cd $h

#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/Original", 1e9, WhichTime="Fix", wantsum=False, wantplot=False, mode="triple", mA=1.105, mB=.934)' >> $d/Original/run.pipe
#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/DiskB-2", 1e7, WhichTime="Fix", wantsum=False, wantplot=False, mode="binary", mA=1.105, mB=.934)' >> $d/DiskB-2/run.pipe
#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/DiskB-3", 1e7, WhichTime="Fix", wantsum=False, wantplot=False, mode="triple", mA=1.105, mB=.934)' >> $d/DiskB-3/run.pipe

	R CMD BATCH -$d '../Analysis/ReadDisk.R'
#	mv ReadDisk.Rout $d
#	python -c 'import Compare;Compare.ComparePipedParams("'$d'",cases=["O","B"])'

#	python -c 'import AlphaCenModule as AC;print(AC.GetLastTime("'$d'/DiskB-2"))'
#	python -c 'import AlphaCenModule as AC;print(AC.GetLastTime("'$d'/DiskB-3"))'

#	echo $d/DiskB-2
#	python -c 'import AlphaCenModule as AC;print(AC.GetAEILastTime("'$d'/DiskB-2"))'
#	cd $d/DiskB-2/Out
#	mkdir -p AeiOutFiles
#	\rm *.aei
#	./elem
#	mv *.aei AeiOutFiles
#	cd $h

#	echo $d/DiskB-3
#	python -c 'import AlphaCenModule as AC;print(AC.GetAEILastTime("'$d'/DiskB-3"))'
#	cd $d/DiskB-3/Out
#	mkdir -p AeiOutFiles
#	\rm *.aei
#	./elem
#	mv *.aei AeiOutFiles
#	cd $h

#	echo '\cp -p '$d/DiskB-2/Out/Backup/* $d/DiskB-2/Out
#	echo '\cp -p '$d/DiskB-3/Out/Backup/* $d/DiskB-3/Out
done

echo 'final i value: '$i
#R CMD BATCH -$i -$hom '../Analysis/CompareProxSurvivals.R'


