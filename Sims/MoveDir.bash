#!/bin/bash
############################################################################### 
### Get name of next available directory
### > ./MoveDir.bash StartDir NewCategory
### NewCategory = 'Err' or 'Proxlike' (where sim should be moved to)

### Determine next available dirname depending on type
name=$(python -c 'import FindNextName; print(FindNextName.NextSubDir("'$2'"))')
echo $name

### Copy current directory to the new place
#\cp -rp $1 $type/$name

