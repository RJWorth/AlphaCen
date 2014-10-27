#!/bin/bash
############################################################################### 
# Copy a proxima-like sim $1 to the Proxlike directory $2
# and create directory for default extension $ext
# e.g: ./ProxSave.bash SDir1

pwd=$(pwd)
echo $pwd

ext='6-3'

### Determine name to copy to
name=$(./MoveDir.bash $1 'Proxlike')
echo 'Copy '$1' to Proxlike/'$name'/Original'

### Copy to Proxlike directory
mkdir Proxlike/$name
\cp -rp $1 Proxlike/$name/Original
mv Proxlike/$name/Original/Out/merc_AC* Proxlike/$name/Original/Out/merc_ext

### Make extension directory
#\cp -rp Proxlike/$name/Original Proxlike/$name/$ext
#\cp -r  Proxlike/param$ext.dmp Proxlike/$2/$ext/Out/param.dmp

