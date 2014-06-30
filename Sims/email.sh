#!/bin/sh
############################################################################### 

machine=$(hostname -s)

echo 'To: rjw274@psu.edu'		>  ACEmail.txt
echo 'From: rjw274@gmail.com'	>> ACEmail.txt
echo 'Subject: re: '$3			>> ACEmail.txt
echo ''							>> ACEmail.txt
echo 'Check '$1' on '$machine', iterations = '$2	>> ACEmail.txt

ssh rjw274@nova.astro.psu.edu 'ssmtp -C ~/SSMTP/default.conf racheljworth@gmail.com < AlphaCen/Sims/ACEmail.txt'
