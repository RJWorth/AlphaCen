#!/bin/bash
############################################################################### 
### Create directories for a new extension
### > ./initialize.bash Proxlike/Prx01 6-3

mkdir gif/$1
mkdir img/$1

if [ ${#2} -gt 0 ]; then
	echo 'subdir'
	mkdir gif/$1/$2
	mkdir img/$1/$2
	mkdir img/$1/$2/gifimgs
else
	echo 'no subdir'
	mkdir img/$1/gifimgs
fi
