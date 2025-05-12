import ROOT
import os,sys
from array import array
import rootUtils as ut
import time
import SndlhcGeo
import re

start_time = time.time()
ROOT.gROOT.SetBatch(ROOT.kTRUE)

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--from_evt", dest="from_evt", required=False, type=int, default=0)
parser.add_argument("--to_evt", dest="to_evt", required=True, type=int, default=None)
parser.add_argument("--numu", dest="numu", required=False, type=bool, default=None)
parser.add_argument("--muDIS", dest="muDIS", required=False, type=bool, default=None)
options = parser.parse_args()

from_evt = options.from_evt
to_evt = options.to_evt

if options.numu:
  pathSim = '/eos/experiment/sndlhc/users/dancc/NUSIM/numu_inBrick21/1'
  pathRec = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_muon1.3E5/b000021'
  pathOut = pathRec+'/eff_study'
elif options.muDIS: ##metti i tuoi path
  pathSim = ''
  pathRec = ''
  pathOut = ''
else:
  print('No input')
  sys.exit(1)

simName = '/sndLHC.Genie-TGeant4.root'
simFile = ROOT.TFile.Open(pathSim+simName)
sTree = simFile.cbmsim

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()

outFileName = pathOut+f'/eff_study.root'
outFile = ROOT.TFile.Open(outFileName, 'RECREATE')
ntuple = ROOT.TNtuple("cbmsim", "Ntuple of nu",'evID:ntracks_true:ttx:tty:nseg_true:vtx_found:ntrks:tracks_found:tx:ty:nseg:s_moth:s_flag:')
ntuple.SetDirectory(outFile)

###### EVENT LOOP ##############
for i_event, event in enumerate(sTree):
  if i_event < from_evt: continue
  if i_event >= to_evt: break
  if i_event%10 == 0: print("Sanity check, current event ", i_event)

  ## opening trk and vtx files per event
  dproc = ROOT.EdbDataProc()
  tproc = ROOT.EdbDataProc()
  gAli = dproc.PVR()
  tAli = dproc.PVR()
  scancond = ROOT.EdbScanCond()
  scancond.SetSigma0(0.5,0.5,0.0025,0.0025)
  scancond.SetDegrad(5)
  gAli.SetScanCond(scancond)
  tAli.SetScanCond(scancond)
  vertexrec = ROOT.EdbVertexRec()
  vertexrec.SetPVRec(gAli)
  vertexrec.eDZmax=3000.
  vertexrec.eProbMin=0.001
  vertexrec.eImpMax=3.5
  vertexrec.eUseMom=False
  vertexrec.eUseSegPar=False
  vertexrec.eQualityMode=0

  trk_cut = 'npl<55'
  vtx_cut = 'flag==0&&vz>-77585&&vz<0'

  trk_file = pathRec+f"/b000021.0.0.{i_event+1}.trk.root"
  vtx_file = pathRec+f"/b000021.0.0.{i_event+1}.vtx.root"
  print("opening file: ", vtx_file)
  if not os.path.isfile(vtx_file):
      print("Vertex file not found, ", vtx_file)
      continue
  tproc.ReadTracksTree(tAli, trk_file, trk_cut)
  dproc.ReadVertexTree(vertexrec, vtx_file, vtx_cut)

  vertices = gAli.eVTX
  tracks = tAli.eTracks

  ## creating graphs for displays
  g1x = ROOT.TGraph()
  g1x.SetName(f"EmulsionDetPoint x {i_event}")
  g2x = ROOT.TGraph()
  g2x.SetName(f"FEDRA tracks x {i_event}")
  g3x = ROOT.TGraph()
  g3x.SetName(f"FEDRA vertex x {i_event}")
  g1y = ROOT.TGraph()
  g1y.SetName(f"EmulsionDetPoint y {i_event}")
  g2y = ROOT.TGraph()
  g2y.SetName(f"FEDRA tracks y {i_event}")
  g3y = ROOT.TGraph()
  g3y.SetName(f"FEDRA vertex y {i_event}")


  ## loop on emudetpoint of charged track

  ## loop on tracks segments

  ## loop on vertex segments
  


  g1x.Write()
  g1y.Write()
  g2x.Write()
  g2y.Write()
  g3x.Write()
  g3y.Write()

  del g1x
  del g2x
  del g3x
  del g1y
  del g2y
  del g3y
  del gAli
  del tAli  

ntuple.Write()
outFile.Write()
outFile.Close() 

print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')