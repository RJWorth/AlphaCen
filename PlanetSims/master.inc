#!/bin/bash
##############################################################################

machine=$(hostname -s)
home=$(pwd)

### Sim parameters
ti=0
tf=3
tDmp=3

### System parameters not changed often
mStar=0.934
rEj=1e2
starA=False
iMax=0.

### Generate a new disk?
new=F	# T or F

### Disk properties
rtr=2.77  #3.08
sigma=3.

### Disk properties, not used yet
rice=2.7
alpha=1.5

### Directory locations (when running code in sim dir)
home=..
c=../Code
in=In/
out=Out/


