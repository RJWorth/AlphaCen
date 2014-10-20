
prf='Err/'
Dirs=(05 07 08 15 24 27)

for i in ${Dirs[*]}
do
	echo '******************* '$i' *******************'

	python -c 'import AlphaCenModule as AC; AC.Summary("'$prf$i'", 1e9, WhichTime="Err", mA=.123,mB=.123, wantsum=False,wantplot=True)'

done
