#!/bin/bash
############################################################################### 
# Copy a proxima-like sim $1 to the Proxlike directory $2
# and create directory for default extension $ext
# e.g: ./ProxSave.bash SDir1 Prx01

pwd=$(pwd)
echo $pwd

ext='6-3'

### Copy to Proxlike directory
mkdir Proxlike/$2
\cp -rp $1 Proxlike/$2/Original
mv Proxlike/$2/Original/Out/merc_AC$1 Proxlike/$2/Original/Out/merc_ext

### Make extension directory
#\cp -rp Proxlike/$2/Original Proxlike/$2/$ext
#\cp -r  Proxlike/param$ext.dmp Proxlike/$2/$ext/Out/param.dmp

