#!/bin/bash
# Sending mail to remote user

sender="rjw274@psu.edu"
receiver="rjw274@psu.edu"
body=$1
subj="AlphaCen notification"

echo $body | mail $receiver -s "$subj"
