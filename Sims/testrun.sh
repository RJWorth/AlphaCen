
#gfortran -w -O1 -o TestA/Out/merc_ACTestA Files/mercury_TG.for
gfortran -w -O1 -o TestA/Out/merc_ACTestA Files/merc_debugging.f

cd TestA/Out
./merc_ACTestA
#./merc_ACSDir1
cd ../..
