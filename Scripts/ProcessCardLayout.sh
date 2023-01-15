#!/bin/sh
#-------------------------------------------------------------------------------
# ProcessCardLayout.sh
#-------------------------------------------------------------------------------
# 2022-10-04
# curtis
#-------------------------------------------------------------------------------
# Turn a CSV describing card layout into a set of files with card layout
#-------------------------------------------------------------------------------

TYPE="$1"
INPUT="$2"

INPUTCARDS="./Cards/FrontBack/$TYPE"
OUTPUT="./Cards/Sheet"

./Scripts/ParseCardLayout.R "$TYPE" "$INPUT" |
    while read l
    do
        ARGLIST=($l)
        ./Scripts/MakeCardFrontBackSheet.sh "$INPUTCARDS/fb_${ARGLIST[1]}.png" \
                                            "$INPUTCARDS/fb_${ARGLIST[2]}.png" \
                                            "$INPUTCARDS/fb_${ARGLIST[3]}.png" \
                                            "$INPUTCARDS/fb_${ARGLIST[4]}.png" > \
        "$OUTPUT/${ARGLIST[0]}sheet.png"
    done

# Mysteriously does not work VV
# img2pdf $(ls $OUTPUT/$TYPE*.png) -o "$OUTPUT/$TYPE.pdf"
