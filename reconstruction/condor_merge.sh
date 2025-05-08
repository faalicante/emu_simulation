#!/bin/bash

export EOSSHIP=root://eospublic.cern.ch/

EVENT=$2
CELL=$3
XPOS=$4
YPOS=$5
OUT_DIR=$6
NU_DIR=/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_eff9_smear1_evt+1
MU_DIR=/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1.3E5/cell_reco/$CELL

echo "Set up SND environment"
SNDBUILD_DIR=/afs/cern.ch/work/f/falicant/public/SNDBUILD/sw
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
echo "Loading FEDRA"
source /afs/cern.ch/work/f/falicant/public/fedra/setup_new.sh

export LD_PRELOAD=/cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/sw/slc9_x86-64/XRootD/latest/lib/libXrdPosixPreload.so
export XROOTD_VMP=eospublic.cern.ch:/eos=/eos

MAIN_DIR=$PWD
cd $MAIN_DIR
MY_DIR=$EVENT
for PLATENUMBER in $(seq 1 60); do
    PLATEFOLDER="$(printf "p%0*d" 3 $PLATENUMBER)"
    mkdir -p ./$MY_DIR/numu/$PLATEFOLDER
    mkdir -p ./$MY_DIR/muon/$PLATEFOLDER
    mkdir -p ./$MY_DIR/b000021/$PLATEFOLDER
    ln -s $NU_DIR/b000021/$PLATEFOLDER/21.$PLATENUMBER.0.0.cp.root ./$MY_DIR/numu/$PLATEFOLDER
    ln -s $MU_DIR/b000021/$PLATEFOLDER/21.$PLATENUMBER.0.0.cp.root ./$MY_DIR/muon/$PLATEFOLDER
done

ln -s $OUT_DIR/neutrino_inbkg.C ./$MY_DIR
mv $MAIN_DIR/track.rootrc ./$MY_DIR/b000021
ln -s $OUT_DIR/vertex.rootrc ./$MY_DIR/b000021

cd $MY_DIR

root -l -q neutrino_inbkg.C\($EVENT,$CELL\)
cd b000021
sed -i "s/XPOS/$XPOS/;s/YPOS/$YPOS/" track.rootrc
makescanset -set=21.0.0.0 -from_plate=60 -to_plate=1 -suff=cp.root -dz=-1315 -v=2 -new
emtra -set=21.0.0.0 -v=2 -new
emvertex -set=21.0.0.0 -v=2

mv b000021.0.0.0.trk.root $MAIN_DIR/b000021.0.0.$(( EVENT + 1)).trk.root
mv b000021.0.0.0.vtx.root $MAIN_DIR/b000021.0.0.$(( EVENT + 1)).vtx.root
