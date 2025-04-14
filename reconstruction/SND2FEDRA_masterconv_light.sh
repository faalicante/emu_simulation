:#!/bin/bash

# Bash script wrapping all of "FairShip2Fedra" conversion tools
# you'll need to copy in this folder the sim file you want to
# convert in FEDRA
SIMFILE=$(ls sndLHC*.root)
GEOFILE=$(ls geofile*.root)
PART=XPART
[ ! -d "$SIMFILE" ] && echo " You are about to convert $SIMFILE "
echo " Bricks and plates folders will be created in the following directory: $(pwd) "
mkdir -p b0000{1..5}{1..4}/p0{01..60}
echo "++++++++++++++++"
echo "FairShip2Fedra parameters are the following:"
cat FairShip2Fedra.rootrc
ulimit -n 1500 #to create many files together
root -l -q fromsndsw2FEDRA.C\(\"$SIMFILE\",\"$GEOFILE\",$PART\)
#echo " Proceeding to doreco.sh "
#source doreco.sh