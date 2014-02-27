#!/bin/sh
# This is some secure program that uses security.

machine=$(hostname -s)
echo $machine

if [ $machine = chloe ]; then
	echo 'chloe match'
elif [ $machine = nova ]; then
	echo 'nova match'
elif [ $machine = myra ]; then
	echo 'myra match'
else
	echo 'no match'
fi


