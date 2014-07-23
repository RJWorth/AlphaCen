#!/bin/sh
############################################################################### 

machine=$(hostname -s)

echo 'To: rjw274@psu.edu'		>  ~/GifEmail.txt
echo 'From: rjw274@gmail.com'	>> ~/GifEmail.txt
echo 'Subject: re: '$2			>> ~/GifEmail.txt
echo ''							>> ~/GifEmail.txt
echo 'Gif '$1' finished.'	>> ~/GifEmail.txt

ssh rjw274@nova.astro.psu.edu 'ssmtp -C ~/SSMTP/default.conf racheljworth@gmail.com < ~/GifEmail.txt'
