#!/bin/bash
###############################################################################
# compare the orbital parameters in the original and disk triple systems

module load python/2.7.3
module load gcc/4.7.1
module load R

h=$(pwd)
echo $h

#pre='Proxlike/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
#Dir=(01 03 05 06 07 09 10 11 12 13 14 15 16 17 19 20 21 22 23 24 26 33)
#Dir=(02 04 08 18 25 27 28 29 30 31 32 34)

#pre='Proxlike/071714/Prx'
#Dir=(1 2 3 4 5 6 7 8 9 10 11 12)


#pre='Proxlike/072314/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26)

pre='Proxlike/081414/Prx'
Dir=(01 02 03 04 05 06 07 08)

for i in ${Dir[*]}; do
	d=$pre$i
	echo '--------------------------------'$d'--------------------------------'
	R CMD BATCH -$d '../Analysis/ReadDisk.R'
	mv ReadDisk.Rout $d
	python -c 'import Compare;Compare.ComparePipedParams("'$d'",cases=["O","B"])'

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

#R CMD BATCH '../Analysis/CompareProxSurvivals.R'


