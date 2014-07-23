#!/bin/bash
############################################################################### 

step1=$(head -1 steps.txt)
steps=$(tail -1 steps.txt)

echo 'first element: '$step1
for x in $steps; 
	do
	echo $x
	done


