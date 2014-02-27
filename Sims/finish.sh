
###############################################################################

# Sum all summary files from all directories

if [ $1 = 1 ]; then
	echo 1
	python -c 'import AlphaCenModule; AlphaCenModule.SumAll(["SDir1","SDir2","SDir3", "CDir1","CDir2","CDir3","CDir4","CDir5", "CDir6","CDir7"],"A")'
fi

if [ $1 = 2 ]; then
	echo 2
	python -c 'import AlphaCenModule; AlphaCenModule.SumAll(["NDir1", "CDir1","CDir2","CDir3","CDir4","CDir5", "CDir6","CDir7", "CDir10", "SDir1","SDir2","SDir3"],"A")'
fi

### Read summary.out into R and make plots
R CMD BATCH --slave plot.R
