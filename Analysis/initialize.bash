#!/bin/bash
############################################################################### 
### Create directories for a new extension
### > ./initialize.bash Proxlike/Prx01 6-3

if [ ${#2} -gt 0 ]; then
	echo 'has subdir'
	mkdir -p gif/$1/$2
#	mkdir img/$1/$2
	mkdir -p img/$1/$2/gifimgs
else
	echo 'no subdir'
	mkdir -p gif/$1
	mkdir -p img/$1/gifimgs
fi
