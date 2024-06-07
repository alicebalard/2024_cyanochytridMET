#!/bin/bash

# Check if exactly two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 title file"
    exit 1
fi

email_subject="$1"
file="$2"

echo "sent from Curta" | mailx -s "$email_subject" -a "$file" "alice.cam.balard@gmail.com" 

