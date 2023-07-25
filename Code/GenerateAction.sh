#!/bin/sh
#-------------------------------------------------------------------------------
# GenerateAction.sh
#-------------------------------------------------------------------------------
# 2022-10-09
# curtis
#-------------------------------------------------------------------------------
# Generate Action cards
#-------------------------------------------------------------------------------

CARDDATA="$1"

exec < "$CARDDATA"
read header
while IFS="," read -r filename title side text
do
    ./Scripts/MakeActionFront.sh "$title" "$side" "$text" > ./Cards/Front/Action/$filename.png
done

