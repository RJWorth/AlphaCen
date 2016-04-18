#!/bin/bash
###############################################################################
# compare the orbital parameters in the original and disk triple systems

#module load python/2.7.3
#module load gcc/4.7.1
#module load R

h=$(pwd)
echo $h

#hom='Proxlike/'
#tag='d'
#pre='Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)

#hom='Proxlike/071714/'
#tag='a'
#pre='Prx'
#Dir=(1 2 3 4 5 6 7 8 9 10 11 12)

#hom='Proxlike/072314/'
#tag='b'
#pre='Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)

hom='Proxlike/081414/'
tag='c'
pre='Prx'
Dir=(01 02 03 04 05 06 07 08)

subdirs=(/Original /DiskA-2 /DiskA-3 /DiskB-2 /DiskB-3)

for i in ${Dir[*]}; do
	d=$hom$pre$i
	echo '--------------------------------'$d'--------------------------------'

	froot=TimeData/TimeData-$tag$i

	if [ ! -f $froot-O.txt ]; then
	scp -rp rjw274@chloe.astro.psu.edu:'/Volumes/Macintosh\ HD\ 2/rjw/AlphaCen/Sims/'$d'/Original/Out/AeiOutFiles/TimeData.txt' TimeData/TimeData-$tag$i-O.txt
	fi

	if [ ! -f $froot-A2.txt ]; then
	scp -rp rjw274@chloe.astro.psu.edu:'/Volumes/Macintosh\ HD\ 2/rjw/AlphaCen/Sims/'$d'/DiskA-2/Out/AeiOutFiles/TimeData.txt' TimeData/TimeData-$tag$i-A2.txt
	fi

	if [ ! -f $froot-A3.txt ]; then
	scp -rp rjw274@chloe.astro.psu.edu:'/Volumes/Macintosh\ HD\ 2/rjw/AlphaCen/Sims/'$d'/DiskA-3/Out/AeiOutFiles/TimeData.txt' TimeData/TimeData-$tag$i-A3.txt
	fi

	if [ ! -f $froot-B2.txt ]; then
	scp -rp rjw274@chloe.astro.psu.edu:'/Volumes/Macintosh\ HD\ 2/rjw/AlphaCen/Sims/'$d'/DiskB-2/Out/AeiOutFiles/TimeData.txt' TimeData/TimeData-$tag$i-B2.txt
	fi

	if [ ! -f $froot-B3.txt ]; then
	scp -rp rjw274@chloe.astro.psu.edu:'/Volumes/Macintosh\ HD\ 2/rjw/AlphaCen/Sims/'$d'/DiskB-3/Out/AeiOutFiles/TimeData.txt' TimeData/TimeData-$tag$i-B3.txt
	fi


#	for j in ${subdirs[*]}; do
#	done

#	cd $d/Original/Out
#	scp -rp rjw274@shapiro.astro.psu.edu:'~/AlphaCen/Sims/'$d'/Original/Out/AeiOutFiles' .
#	scp -p  rjw274@shapiro.astro.psu.edu:'~/AlphaCen/Sims/'$d'/Original/Out/info.out' .
#	./elem
#	mv *.aei AeiOutFiles
#	cd $h

#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/Original", 1e9, WhichTime="Fix", wantsum=False, wantplot=False, mode="triple", mA=1.105, mB=.934)' >> $d/Original/run.pipe
#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/DiskB-2", 1e7, WhichTime="Fix", wantsum=False, wantplot=False, mode="binary", mA=1.105, mB=.934)' >> $d/DiskB-2/run.pipe
#        python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$d'/DiskB-3", 1e7, WhichTime="Fix", wantsum=False, wantplot=False, mode="triple", mA=1.105, mB=.934)' >> $d/DiskB-3/run.pipe

#	R CMD BATCH -$d '../Analysis/ReadDisk.R'
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

#echo 'final i value: '$i
#R CMD BATCH -$i -$hom '../Analysis/CompareProxSurvivals.R'


