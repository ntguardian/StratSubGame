#!/bin/sh
#-------------------------------------------------------------------------------
# GenerateLeverage.sh
#-------------------------------------------------------------------------------
# 2022-09-28
# curtis
#-------------------------------------------------------------------------------
# Generate Leverage cards
#-------------------------------------------------------------------------------

CARDDATA="$1"

exec < "$CARDDATA"
read header
while IFS="," read -r text filename
do
    ./Scripts/MakeLeverageFront.sh "$text" > ./Cards/Front/Leverage/$filename.png
done

