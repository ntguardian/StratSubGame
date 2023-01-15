#!/bin/sh
#-------------------------------------------------------------------------------
# MakeActionFront.sh
#-------------------------------------------------------------------------------
# 2022-09-28
# curtis
#-------------------------------------------------------------------------------
# Make the front of an action card, with plain text
#-------------------------------------------------------------------------------

TITLE=$1
SIDE=$2
TEXT=$3

convert -size 246x378 -background pink xc:pink -fill red -bordercolor Red -border 10x10 -pointsize 18 +size -gravity Center \( -size 226x -background pink pango:"<span size='10240' font='Arial'><span foreground='pink' background='red'><b><i>$TITLE</i></b></span>\n<i>$SIDE</i>: $TEXT</span>" \) -gravity Center -composite png:-

