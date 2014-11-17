#!/bin/bash
###############################################################################
# compare the orbital parameters in the original and disk triple systems

pre='Proxlike/Prx'
#Dir=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34)
Dir=(01 03 05 06 07 09 10 11 12 13 14 15 16 17 19 20 21 22 23 24 26 33)


for i in ${Dir[*]}; do
	d=$pre$i
	echo '--------------------------------'$d'--------------------------------'
	python -c 'import Compare;Compare.ComparePipedParams("'$d'",cases=["O","B"])'
done
