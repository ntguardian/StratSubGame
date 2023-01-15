#!/bin/zsh
#-------------------------------------------------------------------------------
# GenerateCommand.zsh
#-------------------------------------------------------------------------------
# 2022-10-16
# curtis
#-------------------------------------------------------------------------------
# Generate a set of command cards
#-------------------------------------------------------------------------------

CARDDATA="$1"

exec < "$CARDDATA"
read header
while IFS="," read -r mission subtype filename
do
    convert =(./Scripts/MakeCommand.sh "$mission") =(./Scripts/MakeCommand.sh "$subtype") +append -bordercolor White -border 2x2 "./Cards/FrontBack/Command/fb_$filename.png"
done

