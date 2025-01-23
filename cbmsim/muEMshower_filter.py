import os
# Code snippet from Simona. This should go somewhere where it can be used by several different pieces of code (monitoring, analysis, etc)
import pickle

import numpy as np
import ROOT
import rootUtils as ut
import SndlhcGeo
from array import array

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", dest = "inputFile", required = False)
parser.add_argument("-o", dest = "outDir", required = False, default='./')
parser.add_argument("-g", dest = "geoFile", required = False)
parser.add_argument("--filter", dest="filter", required=False, default=False, action='store_true')
parser.add_argument("--ana", dest="ana", required=False, default=False, action='store_true')
args = parser.parse_args()

geoFile = '/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_z289.374023_BRICK11/9790422/geofile_full.Ntuple-TGeant4.root'

# Set up TTrees
treeName = "cbmsim"
ch = ROOT.TChain(treeName)
ch.Add(args.inputFile)

if ch.GetEntries() == 0 :
    print("Chain is empty. Exitting")
    exit(-1)

if args.geoFile:
  geoFile = args.geoFile

snd_geo = SndlhcGeo.GeoInterface(geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

isMC = True
obj = {}
h = {}

def filterFile():
  ch.GetEntry(0)
  f = ch.GetFile()

  tmp = args.inputFile.split('/')
  infilename = tmp[len(tmp)-1]
  tag = '_EMshower'
  output_file = ROOT.TFile(args.outDir+'/'+infilename[:-5]+tag+'.root', "RECREATE")
  output_tree = ch.CloneTree(0)
  # Copy branch list
  branch_list = f.Get("BranchList")
  branch_list_copy = branch_list.Clone()
  branch_list_copy.Write("BranchList", 1)

  HAD_procIDs = [23,13,25,24,26, 27,6]
  interest_procIDs = [5, 7, 8, 9]
  passed_evts = 0

  for event in ch:
    HAD_event = False
    accept_event = False
    for track in event.MCTrack:
      if track.GetProcID() in HAD_procIDs:
        HAD_event = True
        break
      if track.GetEnergy() < 5.: continue
      if track.GetMotherId()==0 and track.GetProcID() in interest_procIDs:
        accept_event = True

    if HAD_event: continue
    if accept_event:
      output_tree.Fill()
      passed_evts +=1

  print(passed_evts, 'have passed the selection')
  output_file.Close()
  output_file.Write()

def analyseEM():
  hists = {}
  ### HISTOS HERE
  ut.bookHist(hists, "InteractionPlaneSF", "Scifiplane with most e^{-}; Scifiplane", 5, 0.5, 5.5)
  ut.bookHist(hists, "InteractionPlaneSF_prim", "Scifiplane with most primary e^{-}; Scifiplane", 5, 0.5, 5.5)
  ut.bookHist(hists, "InteractionPlaneEM", "Wall with most e^{-}; Wall", 5, 0.5, 5.5)
  ut.bookHist(hists, "Nsf_hits", "N of Scifi hits", 600, 0, 600)
  ut.bookHist(hists, "muEn", "Muon Energy", 480, 0, 4800)
  ut.bookHist(hists, "MaxElWall", "Wall vs n. max e^{-}", 5, 0.5, 5.5, 50, 0, 500)
  ut.bookHist(hists, "nEl_muEn", "N max e^{-} vs E_{#mu}", 50, 0, 500, 480, 0, 4800)
  ut.bookHist(hists, "nsf_muEn", "N of Scifi hits vs E_{#mu}", 600, 0, 600, 480, 0, 4800)
  ut.bookHist(hists, "Nsf_afterInt", "N of Scifi hits after the interaction wall;Nscifi", 600, 0, 600)
  ut.bookHist(hists, "Nsf_Nel_afterInt", "N of Scifi hits vs N. e^{-} after the interaction wall;Nscifi;N e^{-}", 600, 0, 600, 50, 0, 500)

  ut.bookHist(hists, "STD_showerX", "STD of shower in X direction", 500, 0, 50)
  ut.bookHist(hists, "STD_showerY", "STD of shower in Y direction", 500, 0, 50)
  ut.bookHist(hists, "STD_showerXY", "STD of shower in XY direction", 500, 0, 50, 500, 0, 50)

  for plane in range(1, 6):
    for inplane in range(1, 6):
      if inplane <= plane:
        ut.bookHist(hists, 'Nsfplane'+str(plane)+'_after_Int_'+str(inplane), "N of Scifi hits in plane "+str(plane)+" after Int in"+str(inplane), 600, 0, 600)

  ### LOOP HERE
  vLeft, vRight = ROOT.TVector3(), ROOT.TVector3()
  fitted_meansX_list = list()
  fitted_meansY_list = list()
  fitted_stdX_list = list()
  fitted_stdY_list = list()

  for ievent, event in enumerate(ch):
    Wall_el = {1:[], 2: [], 3:[], 4:[], 5:[]}
    SF_el = {1:[], 2: [], 3:[], 4:[], 5:[]}
    SF_el_prim = {1:[], 2: [], 3:[], 4:[], 5:[]}
    wall_nel_pair = []
    SF_nel_pair = []
    SF_nel_pair_prim = []
    """for point in event.EmulsionDetPoint:
      if int(point.GetDetectorID()/1e4)==1 and point.GetTrackID()==0:
        muE = ROOT.TMath.Sqrt(point.GetPx()*point.GetPx()+point.GetPy()*point.GetPy()+point.GetPz()*point.GetPz())
        break

    for point in event.EmulsionDetPoint:
      pdg = point.PdgCode()
      trackID = point.GetTrackID()
      Wall = int(point.GetDetectorID()/1e4)
      #print(pdg, trackID, Wall)
      if trackID < 0: continue
      if abs(pdg) != 11: continue
      if trackID in Wall_el[Wall]: continue
      Wall_el[Wall].append(trackID)
    for wall in Wall_el.keys():
      wall_nel_pair.append([wall, len(Wall_el[wall])])
    maxpairEM = max(wall_nel_pair, key=lambda item:item[1])"""
    muE = event.MCTrack[0].GetEnergy()
    hists['muEn'].Fill(muE)
    for point in event.ScifiPoint:
      pdg = point.PdgCode()
      trackID = point.GetTrackID()
      station = point.station()
      if trackID < 0: continue
      if abs(pdg) != 11: continue
      if trackID in Wall_el[station]: continue
      SF_el[station].append(trackID)
      if event.MCTrack[trackID].GetMotherId()==0:
        SF_el_prim[station].append(trackID)
    for wall in SF_el.keys():
      SF_nel_pair.append([wall, len(SF_el[wall])])
    maxpair = max(SF_nel_pair, key=lambda item:item[1])
    if maxpair[1] == 0: continue
    for plane in SF_el_prim.keys():
      SF_nel_pair_prim.append([plane, len(SF_el_prim[plane])])
    maxpair_prim = max(SF_nel_pair_prim, key=lambda item:item[1])
    if len(event.ScifiPoint) == 0: continue
    #print('event', ievent, 'Max of electrons in', maxpair[0], 'n electrons', maxpair[1])
    hists['InteractionPlaneSF'].Fill(maxpair[0])
    hists['InteractionPlaneSF_prim'].Fill(maxpair_prim[0])
    #hists['InteractionPlaneEM'].Fill(maxpairEM[0])
    hists['MaxElWall'].Fill(maxpair[0], maxpair[1])
    hists['nEl_muEn'].Fill(maxpair[1], muE)

    Nsf     = 0
    nsf_statID = {1:0, 2:0, 3:0, 4:0, 5:0}
    Scifi_HitCollectionX = {1:[], 2:[], 3:[], 4:[], 5:[]}
    Scifi_HitCollectionY = {1:[], 2:[], 3:[], 4:[], 5:[]}
    for point in event.ScifiPoint:
      station = point.station()
      if station != maxpair[0]: continue
      if abs(point.PdgCode()) != 13: continue
      if point.GetTrackID()!=0: continue
      track_conv_pos = ROOT.TVector3(point.GetX(), point.GetY(), point.GetZ())
      break

    for aHit in event.Digi_ScifiHits:
      if not aHit.isValid(): continue
      station = aHit.GetStation()
      detID = aHit.GetDetectorID()
      scifiDet.GetSiPMPosition(detID, vLeft, vRight)
      Nsf+=1
      nsf_statID[station] +=1
      if station == maxpair[0]:
        if aHit.isVertical():
          Scifi_HitCollectionX[station].append(vLeft.X())
        else:
          Scifi_HitCollectionY[station].append(vRight.Y())

    xpositions = np.array([pos for pos in Scifi_HitCollectionX[maxpair[0]] if pos > track_conv_pos.X()-10 and pos < track_conv_pos.X()+10])
    ypositions = np.array([pos for pos in Scifi_HitCollectionY[maxpair[0]] if pos > track_conv_pos.Y()-10 and pos < track_conv_pos.Y()+10])

    if len(xpositions) == 0 or len(ypositions) == 0: continue
    fitted_meansX_list.append(np.mean(xpositions))
    fitted_meansY_list.append(np.mean(ypositions))
    fitted_stdX_list.append(np.std(xpositions))
    fitted_stdY_list.append(np.std(ypositions))


    hists['STD_showerX'].Fill(np.std(xpositions))
    hists['STD_showerY'].Fill(np.std(ypositions))
    hists['STD_showerXY'].Fill(np.std(xpositions), np.std(ypositions))

    hists['Nsf_hits'].Fill(Nsf)
    hists['nsf_muEn'].Fill(Nsf, muE)
    hists['Nsf_afterInt'].Fill(nsf_statID[maxpair[0]])
    hists['Nsf_Nel_afterInt'].Fill(nsf_statID[maxpair[0]], maxpair[1])
    for plane in range(1, 6):
      if plane >= maxpair[0]:
        hists['Nsfplane'+str(plane)+'_after_Int_'+str(maxpair[0])].Fill(nsf_statID[plane])

  fout = ROOT.TFile.Open(args.outDir+'/EMhists_Sept24.root', "RECREATE")
  for _h in hists.values():
    _h.Write()
  fout.Write()
  fout.Close()


if args.filter: filterFile()
if args.ana: analyseEM()