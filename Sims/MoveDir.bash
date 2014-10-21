#!/bin/bash
############################################################################### 
### Rename and move Prox or Err directories as appropriate
### > ./MoveDir.bash StartDir NewDir
### NewDir = 'Err' or 'Proxlike' (where sim should be moved to)

### Determine next available dirname depending on type
name=$(python -c 'import FindNextName; print(FindNextName.NextSubDir("'$2'"))')
echo $name

### Copy current directory to the new place
#\cp -rp $1 $type/$name

