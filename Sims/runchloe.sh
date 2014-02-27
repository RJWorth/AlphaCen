
############################################################################### 
# Start time
t1=$(date +%s)

### clean?

# will maybe need a separate starter script for each dir, since they loop
# (or one that starts up several separate scripts that run in parallel) 
## Loop over directories

cent=A

aBmin=23.4
aBmax=28.1
eBmin=0.0
eBmax=0.52
iBmin=0.0
iBmax=0.0

aCmin=1.0	# *aB(1+eB)/(1-eC), i.e. min(periapse of C) = this*(apoapse of B)
aCmax=5.0	# ditto, but max(periapse)
eCmin=0.0
eCmax=0.99
iCmin=0.0
iCmax=180

echo ../In/big.in         >  $1/Out/files.in
echo ../In/small.in       >> $1/Out/files.in
echo ../../Files/param.in >> $1/Out/files.in
#echo ../In/param.in      >> $1/Out/files.in
echo xv.out               >> $1/Out/files.in
echo ce.out               >> $1/Out/files.in
echo info.out             >> $1/Out/files.in
echo runlog.out           >> $1/Out/files.in
echo big.dmp              >> $1/Out/files.in
echo small.dmp            >> $1/Out/files.in
echo param.dmp            >> $1/Out/files.in
echo restart.dmp          >> $1/Out/files.in



### Repeat simulation n times and record summary of collisions
for j in {1..100}
#for j in $(seq 1 100)
do
	echo 'run: '$1', '$j
	\rm $1/Out/*.dmp
	\rm $1/Out/*.out
#### Randomize B and C parameters (a, e, i)
	python -c 'import AlphaCenModule; AlphaCenModule.MakeBigRand( "'$1'",'$j',"'$cent'", '$aBmin','$aBmax','$eBmin','$eBmax','$iBmin','$iBmax', '$aCmin','$aCmax','$eCmin','$eCmax','$iCmin','$iCmax')'

### Randomize rock positions
# currently no rocks needed

#### Run mercury
	gfortran -w -o $1/Out/merc_debugging Files/merc_debugging.f
	cd $1/Out;	./merc_debugging;	cd ../..

### If needing element:
	gfortran-4.2 -w -o $1/Out/elem Files/elem.for
	cd $1/Out;	./elem;	cd ../..
	\mv $1/Out/*.aei $1/Out/AeiOutFiles

### Write summary
	python -c 'import AlphaCenModule; AlphaCenModule.Summary("'$1'","'$j'","'$cent'")'


done

# Stop time:
t2=$(date +%s)
echo $1"	$(echo "$t2 - $t1"|bc )" >> runtime.txt

