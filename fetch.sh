#!/bin/bash

FILES=$1
USER=$2
PASS=$3

for LINE in `cat $FILES`
do
    URL=$(echo $LINE | cut -f1 -d",")
    DEST=$(echo $LINE | cut -f2 -d",")
    
    echo "====================================================================="
    echo "Fetching BAM file for $DEST"
    echo "====================================================================="
    if [ -f "$DEST" ]
    then
        echo "---------------------------------------------------------------------"
        echo "$DEST exists! Skipping..."
        echo "---------------------------------------------------------------------"
    else
        curl --user $USER":"$PASS \
            --header "Content-Type: application/bam" \
            --location $URL > $DEST
    fi
    echo " "
done


