#!/bin/sh
#-------------------------------------------------------------------------------
# MakeIntelFront.sh
#-------------------------------------------------------------------------------
# 2022-09-28
# curtis
#-------------------------------------------------------------------------------
# Make the front of an intel card, with plain text
#-------------------------------------------------------------------------------

TEXT=$1

convert -background lightgreen -fill green -bordercolor green -border 10x10 -pointsize 18 -gravity Center -size 241x373 caption:"$TEXT" png:-

