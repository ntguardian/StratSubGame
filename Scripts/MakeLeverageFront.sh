#!/bin/sh
#-------------------------------------------------------------------------------
# MakeLeverageFront.sh
#-------------------------------------------------------------------------------
# 2022-09-28
# curtis
#-------------------------------------------------------------------------------
# Make the front of a leverage card, with plain text
#-------------------------------------------------------------------------------

TEXT=$1

convert -background lightblue -fill blue -bordercolor Blue -border 10x10 -pointsize 18 -gravity Center -size 246x378 caption:"$TEXT" png:-

