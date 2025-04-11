:#!/bin/bash

# Bash script wrapping all of "FairShip2Fedra" conversion tools
# you'll need to copy in this folder the sim file you want to
# convert in FEDRA
RED=$'\033[0;31m' #Red coloring                                                                                                                                                  
NC=$'\033[0m' #No coloring

SIMFILE=/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_BRICK11_RUN1/8461734
GEOFILE=$(ls geofile*.root)
#[ ! -d "$SIMFILE" ] && echo "${RED}++${NC} You are about to convert $SIMFILE ${RED}++${NC}"
echo "${RED}++${NC} Bricks and plates folders will be created in the following directory: $(pwd) ${RED}++${NC}"
#cp /afs/cern.ch/work/f/falicant/public/macros-snd/FEDRA/fromsndsw2FEDRA.C .
for i in $(seq 1 5)
	do
		mkdir b0000${i}{1..4}
		mkdir b0000${i}{1..4}/p0{01..60}
	done
echo "++++++++++++++++"
echo "FairShip2Fedra parameters are the following:"
cat FairShip2Fedra.rootrc
ulimit -n 1500 #to create many files together
root -l -q fromsndsw2FEDRA.C\(\"$SIMFILE\",\"$GEOFILE\",02)1)\)
echo "${RED}++${NC} Proceeding to doreco.sh ${RED}++${NC}"
#source doreco.sh
