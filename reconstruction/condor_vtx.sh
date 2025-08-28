#!/bin/bash

export EOSSHIP=root://eospublic.cern.ch/

EVENT=$2
OUT_DIR=$3

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

mkdir -p ./$MY_DIR/b000021
ln -s $OUT_DIR/b000021/b000021.0.0.0.set.root ./$MY_DIR/b000021/
ln -s $OUT_DIR/b000021/b000021.0.0.$(( EVENT + 1)).trk.root ./$MY_DIR/b000021/b000021.0.0.0.trk.root
ln -s $OUT_DIR/vertex_disc.rootrc ./$MY_DIR/b000021
ln -s $OUT_DIR/vertex_edi.rootrc ./$MY_DIR/b000021

cd $MY_DIR/b000021

cp vertex_disc.rootrc vertex.rootrc
emvertex -set=21.0.0.0 -v=1
emvertex -set=21.0.0.0 -v=1 -r

cp b000021.0.0.0.vtx.root $MAIN_DIR/b000021.0.0.$(( EVENT + 1)).vtx.root
cp b000021.0.0.0.vtx.discimp.root $MAIN_DIR/b000021.0.0.$(( EVENT + 1)).vtx.discimp.root

ln -s b000021.0.0.0.vtx.discimp.root b000021.0.0.0.vtx.root
cp vertex_edi.rootrc vertex.rootrc
emvertex -set=21.0.0.0 -v=1 -r -fit

cp b000021.0.0.0.vtx.refit.root $MAIN_DIR/b000021.0.0.$(( EVENT + 1)).vtx.refit.root
