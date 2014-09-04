
###############################################################################

# Sum all summary files from all directories

#if [ $1 = 1 ]; then
#	echo 1
	python -c 'import AlphaCenModule; AlphaCenModule.SumAll(["SDir1","SDir2","SDir3","SDir4","SDir5", "SDir6","SDir7", "CDir1","CDir2","CDir3","CDir4","CDir5","CDir6","CDir7","CDir8","CDir9","CDir10" ],"A")'
#fi

### Read summary.out into R and make plots
R CMD BATCH -1 plot.R Plots/plot.Rout

