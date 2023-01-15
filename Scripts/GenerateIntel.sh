#!/bin/sh
#-------------------------------------------------------------------------------
# GenerateIntel.sh
#-------------------------------------------------------------------------------
# 2022-09-29
# curtis
#-------------------------------------------------------------------------------
# Generate a set of intel cards from input data
#-------------------------------------------------------------------------------

CARDDATA="$1"

exec < "$CARDDATA"
read header
while IFS="," read -r text filename
do
    ./Scripts/MakeIntelFront.sh "$text" > ./Cards/Front/Intel/$filename.png
done

