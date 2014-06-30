#!/bin/bash
############################################################################### 

echo $2
end=${2: -1}
echo $end

echo ../In/big.in         >  $1/Out/files.in
echo ../In/small.in       >> $1/Out/files.in
#echo ../../Files/param.in >> $1/Out/files.in	#one master param.in
echo ../In/param.in       >> $1/Out/files.in	#individual param.in
echo xv.out               >> $1/Out/files.in
echo ce.out               >> $1/Out/files.in
echo info.out             >> $1/Out/files.in
	if [ $end = 'f' ]; then
	echo runlog.out       >> $1/Out/files.in
	fi
echo big.dmp              >> $1/Out/files.in
echo small.dmp            >> $1/Out/files.in
echo param.dmp            >> $1/Out/files.in
echo restart.dmp          >> $1/Out/files.in

