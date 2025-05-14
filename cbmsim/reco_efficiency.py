import ROOT
import fedrarootlogon
import os,sys
from array import array
import time
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

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

from_evt = options.from_evt
to_evt = options.to_evt

if options.numu:
  pathSim = '/eos/experiment/sndlhc/users/dancc/NUSIM/numu_inBrick21/1'
  # pathRec = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_muon1.3E5/b000021'
  pathRec = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numu_inBrick21/eff100_nosmear/b000021'
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

colors = [1,2,4,6,7,8,9,13,28,40,42,46,38,30]

outFileName = pathOut+f'/eff_study.root'
outFile = ROOT.TFile.Open(outFileName, 'RECREATE')
ntuple = ROOT.TNtuple("cbmsim", "Ntuple of nu",'evID:ntracks_true:ttx:tty:nseg_true:vtx_found:ntrks:ff:ip:prob:tracks_found:tx:ty:nseg')
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
  vtx_cut = 'vz>-77585&&vz<0'

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
  g1x.SetName(f"EmuDetPointx_{i_event}")
  g1y = ROOT.TGraph()
  g1y.SetName(f"EmuDetPointy_{i_event}")
  g2xl = []
  g2yl = []
  g4xl = []
  g4yl = []

  ## loop on emudetpoint of charged track
  emuTracks = {}
  for p in event.EmulsionDetPoint:
    trackID = p.GetTrackID()
    if trackID < 0: continue
    if not MyPDG.GetParticle(p.PdgCode()): continue
    if MyPDG.GetParticle(p.PdgCode()).Charge()==0: continue
    if int(p.GetDetectorID()/1000) != 21: continue
    if event.MCTrack[trackID].GetMotherId()==0:
      AddToDict(emuTracks, trackID)
      ex = p.GetX() * 1E+4 + 473000
      ey = p.GetY() * 1E+4 - 158000
      ez = p.GetZ()
      g1x.SetPoint(g1x.GetN(), ez, ex)
      g1y.SetPoint(g1y.GetN(), ez, ey)
  ntracks_true = len(emuTracks)

  ## loop on tracks segments
  sxlistsig = []
  sylistsig = []
  szlistsig = []
  for track in tracks:
    nseg = track.N()
    sxlist = []
    sylist = []
    szlist = []
    draw = False
    for iseg in range(nseg):
      seg = track.GetSegment(iseg)
      signal = seg.Flag()
      mother = seg.Aid(0)
      sxlist.append(seg.X())
      sylist.append(seg.Y())
      szlist.append(seg.Z())
      if signal == 1 and mother == 0:
        sxlistsig.append(seg.X())
        sylistsig.append(seg.Y())
        szlistsig.append(seg.Z())
        draw = True
    if draw:
      sxArr = array('d', sxlist)
      syArr = array('d', sylist)
      szArr = array('d', szlist)
      g2xl.append(ROOT.TGraph(len(sxlist), szArr, sxArr))
      # g2xl[len(g2xl)-1].SetName(f"FEDRAtx_{len(g2xl)-1}_{i_event}")
      g2yl.append(ROOT.TGraph(len(sylist), szArr, syArr))
      # g2yl[len(g2yl)-1].SetName(f"FEDRAty_{len(g2yl)-1}_{i_event}")
  sxlistsigArr = array('d', sxlistsig)
  sylistsigArr = array('d', sylistsig)
  szlistsigArr = array('d', szlistsig)
  g3x = ROOT.TGraph(len(sxlistsig), szlistsigArr, sxlistsigArr)
  g3x.SetName(f"FEDRAtx_{i_event}")
  g3y = ROOT.TGraph(len(sylistsig), szlistsigArr, sylistsigArr)
  g3y.SetName(f"FEDRAty_{i_event}")

  ## loop on vertex segments
  for vertex in vertices:
    n = vertex.N()
    sxlist = []
    sylist = []
    szlist = []
    draw = False
    for itrack in range(n):
      track = vertex.GetTrack(itrack)
      nseg = track.N()
      for iseg in range(nseg):
        seg = track.GetSegment(iseg)
        signal = seg.Flag()
        mother = seg.Aid(0)
        sxlist.append(seg.X())
        sylist.append(seg.Y())
        szlist.append(seg.Z())
        if signal == 1 and mother == 0:
          draw = True
    if draw:
      print("found vertex")
      sxArr = array('d', sxlist)
      syArr = array('d', sylist)
      szArr = array('d', szlist)
      g4xl.append(ROOT.TGraph(len(sxlist), szArr, sxArr))
      g4yl.append(ROOT.TGraph(len(sylist), szArr, syArr))
  found_vtx = len(g4xl)

  ## drawing graphs
  c1 = ROOT.TCanvas(f"c1_{i_event}", f"c1_{i_event}", 1500, 1000)
  c2 = ROOT.TCanvas(f"c2_{i_event}", f"c2_{i_event}", 1500, 1000)
  c1.Divide(2, 2)
  c2.Divide(2, 2)

  c1.cd(1)
  g1x.SetTitle(f"Event {i_event} - EmulsionDetPoint x;z;x")
  g1x.SetMarkerStyle(20)
  g1x.SetMarkerColor(ROOT.kBlack)
  g1x.SetMarkerSize(1)
  gxMax = g1x.GetYaxis().GetXmax()
  gxMin = g1x.GetYaxis().GetXmin()
  g1x.GetYaxis().SetRangeUser(gxMin, gxMax)
  print(gxMin, gxMax)
  g1x.Draw("AP")

  c1.cd(3)
  g1y.SetTitle(f"Event {i_event} - EmulsionDetPoint y;z;y")
  g1y.SetMarkerStyle(20)
  g1y.SetMarkerColor(ROOT.kBlack)
  g1y.SetMarkerSize(1)
  gyMax = g1y.GetYaxis().GetXmax()
  gyMin = g1y.GetYaxis().GetXmin()
  g1y.GetYaxis().SetRangeUser(gyMin, gyMax)
  g1y.Draw("AP")

  g2x = ROOT.TMultiGraph()
  g2x.SetTitle(f"Event {i_event} - FEDRA tracks x;z;x")
  for ig, gx in enumerate(g2xl):
    gx.SetMarkerStyle(20)
    gx.SetMarkerColor(colors[ig%len(colors)])
    gx.SetMarkerSize(1)
    g2x.Add(gx)
  c1.cd(2)
  g2x.GetYaxis().SetRangeUser(gxMin, gxMax)
  g2x.Draw("AP")
  c2.cd(1)
  g2x.Draw("AP")

  g2y = ROOT.TMultiGraph()
  g2y.SetTitle(f"Event {i_event} - FEDRA tracks y;z;y")
  for ig, gy in enumerate(g2yl):
    gy.SetMarkerStyle(20)
    gy.SetMarkerColor(colors[ig%len(colors)])
    gy.SetMarkerSize(1)
    g2y.Add(gy)
  c1.cd(4)
  g2y.GetYaxis().SetRangeUser(gyMin, gyMax)
  g2y.Draw("AP")
  c2.cd(3)
  g2y.Draw("AP")

  c2.cd(2)
  g3x.SetTitle(f"Event {i_event} - FEDRA sig tracks x;z;x")
  g3x.SetMarkerStyle(20)
  g3x.SetMarkerColor(ROOT.kRed)
  g3x.SetMarkerSize(1)
  g3x.GetYaxis().SetRangeUser(gxMin, gxMax)
  g3x.Draw("AP")

  c2.cd(4)
  g3y.SetTitle(f"Event {i_event} - FEDRA sig tracks y;z;y")
  g3y.SetMarkerStyle(20)
  g3y.SetMarkerColor(ROOT.kRed)
  g3y.SetMarkerSize(1)
  g3y.GetYaxis().SetRangeUser(gyMin, gyMax)
  g3y.Draw("AP")

  c1.Update()
  c2.Update()

  outFile.cd()
  c1.SaveAs(pathOut+f"/c1_{i_event}.png")
  c2.SaveAs(pathOut+f"/c2_{i_event}.png")
  c1.Write()
  c2.Write()

  if found_vtx > 0:
    c3 = ROOT.TCanvas(f"c3_{i_event}", f"c3_{i_event}", 1500, 1000)
    c3.Divide(2, 2)

    c3.cd(2)
    g3x.Draw("AP")
    c3.cd(4)
    g3y.Draw("AP")
    
    c3.cd(1)
    g4x = ROOT.TMultiGraph()
    g4x.SetTitle(f"Event {i_event} - FEDRA vertex x;z;x")
    for ig, gx in enumerate(g4xl):
      gx.SetMarkerStyle(20)
      gx.SetMarkerColor(colors[ig%len(colors)])
      gx.SetMarkerSize(1)
      g4x.Add(gx)
    g4x.GetYaxis().SetRangeUser(gxMin, gxMax)
    g4x.Draw("AP")
  
    c3.cd(3)
    g4y = ROOT.TMultiGraph()
    g4y.SetTitle(f"Event {i_event} - FEDRA vertex x;z;y")
    for ig, gy in enumerate(g4yl):
      gy.SetMarkerStyle(20)
      gy.SetMarkerColor(colors[ig%len(colors)])
      gy.SetMarkerSize(1)
      g4y.Add(gy)
    g4y.GetYaxis().SetRangeUser(gyMin, gyMax)
    g4y.Draw("AP")
    c3.Update()
    outFile.cd()
    c3.Write()
    c3.SaveAs(pathOut+f"/c3_{i_event}.png")

  del gAli
  del tAli

outFile.cd()
ntuple.Write()
outFile.Write()
outFile.Close() 

print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')