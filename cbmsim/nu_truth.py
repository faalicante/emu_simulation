import ROOT
import os,sys
from array import array
import rootUtils as ut
import time
from datetime import date
today = date.today().strftime('%d%m%y')
start_time = time.time()
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def getOriginAndDims(nodepath):
  nav = ROOT.gGeoManager.GetCurrentNavigator()
  nav.cd(nodepath)
  O = {'X':0, 'Y':0, 'Z':0}
  D = {'X':0, 'Y':0, 'Z':0}
  N = nav.GetCurrentNode()
  S = N.GetVolume().GetShape()
  D['X'], D['Y'], D['Z'] = S.GetDX(),S.GetDY(),S.GetDZ()
  O['X'], O['Y'], O['Z'] = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
  OriginArray = array('d', O.values())
  OriginTrans = array('d', [0, 0, 0])
  nav.LocalToMaster(OriginArray, OriginTrans)
  O['X'], O['Y'], O['Z'] = OriginTrans[0],OriginTrans[1],OriginTrans[2]
  return O, D

def getWallRanges():
  Walls = {}
  Ranges = {'X': list(), 'Y': list(), 'Z': list()}
  for i in range(5):
    node = '/cave_1/Detector_0/volTarget_1/Wall_{0}'.format(i)
    Walls[i] = getOriginAndDims(node)
    for proj in ['X', 'Y', 'Z']:
      Ranges[proj].append([Walls[i][0][proj]-Walls[i][1][proj], Walls[i][0][proj]+Walls[i][1][proj]])
  return Ranges

def getBrickRanges():
  Bricks = {}
  Ranges = {'X': list(), 'Y': list(), 'Z': list()}
  for wall in range(5):
    Bricks[wall] = {}
    for row in range(2):
      Bricks[wall][row] = {}
      for brick in range(1,-1,-1):
        node = '/cave_1/Detector_0/volTarget_1/Wall_{0}/Row_{1}/Brick_{2}'.format(wall, row, brick)
        Bricks[wall][row][brick] = getOriginAndDims(node)
        for proj in ['X', 'Y', 'Z']:
          Ranges[proj].append([Bricks[wall][row][brick][0][proj]-Bricks[wall][row][brick][1][proj], Bricks[wall][row][brick][0][proj]+Bricks[wall][row][brick][1][proj]])
  return Ranges

def isInWall(pos, wall, ranges):
  if pos.X() > ranges['X'][wall][0] and pos.X() < ranges['X'][wall][1]:
    if pos.Y() > ranges['Y'][wall][0] and pos.Y() < ranges['Y'][wall][1]:
      if pos.Z() > ranges['Z'][wall][0] and pos.Z() < ranges['Z'][wall][1]:
        return True

def getBrickInt(vtx, brickRanges):
  ret = False
  for brick_int in range(20):
    if isInWall(vtx, brick_int, brickRanges):
      ret = True
      break
  return ret, brick_int

def getWallInt(vtx, wallRanges):
  ret = False
  for wall_int in range(5):
    if isInWall(vtx, wall_int, wallRanges):
      ret = True
      break
  return ret, wall_int+1

def decodeBrick(brick):
  wall = brick // 4 + 1
  brick = wall*10 + (brick % 4 + 1)
  return wall, brick

def FindMuons(MCTrack):
    muons = {'In': None, 'Out': None}
    muons['Out'] = 0
    for itrk, track in enumerate(MCTrack):
        if ROOT.TMath.Abs(track.GetPdgCode())!=13: continue
        if track.GetMotherId()!=0:continue
        if track.GetPz() > 0: continue
        muons['In'] = itrk
    return muons


pathSim = '/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/nu_sim_activeemu_withcrisfiles_25_July_2022/'
geoFile =  pathSim + '/geofile_full.Genie-TGeant4.root'
outPath = '/afs/cern.ch/work/f/falicant/public/emu_simulation/out'

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(geoFile)

simName = "/sndLHC.Genie-TGeant4.root"
simFile = ROOT.TFile.Open(pathSim+simName)
sTree = simFile.cbmsim
gst = simFile.gst
nentries = sTree.GetEntries()
print('There are ', nentries, 'entries')

up_to_ev = None
if len(sys.argv)>1:
  up_to_ev = int(sys.argv[1])
  if len(sys.argv) == 3:
    from_ev = int(sys.argv[1])
    to_ev = int(sys.argv[2])
    print('Starting from ev', from_ev, 'to ev', to_ev)

wallRanges = getWallRanges()
brickRanges = getBrickRanges()

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()

h               = {}
e_cuts          = {200:'_e200', 500:'_e500', 1000:'_e1000'}
angle_cuts      = {0:'', 1:'_angle1'}
vis_cuts        = {**angle_cuts, **e_cuts}
for vcut in vis_cuts.values():
  for nu in ['numu', 'nue', 'numu_c', 'nue_c', 'neu']:
    ut.bookHist(h, f'{nu}_energy{vcut}', f';E_{nu}', 240, 0, 4800)
    ut.bookHist(h, f'{nu}_TX{vcut}', ';TX', 300, -1.5, 1.5)
    ut.bookHist(h, f'{nu}_TY{vcut}', ';TY', 300, -1.5, 1.5)
    ut.bookHist(h, f'{nu}_TXTY{vcut}', ';TX;TY', 300, -1.5, 1.5, 300, -1.5, 1.5)
    ut.bookHist(h, f'{nu}_vx{vcut}', f";{nu} vtx X[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
    ut.bookHist(h, f'{nu}_vy{vcut}', f";{nu} vtx Y[cm]", int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
    ut.bookHist(h, f'{nu}_vz{vcut}', f";{nu} vtx Z[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20)
    ut.bookHist(h, f'{nu}_vxy{vcut}', f";{nu} vtx X[cm];{nu} vtx Y[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
    ut.bookHist(h, f'{nu}_vxz{vcut}', f";{nu} vtx Z[cm];{nu} vtx X[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
    ut.bookHist(h, f'{nu}_vyz{vcut}', f";{nu} vtx Z[cm];{nu} vtx Y[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
    ut.bookHist(h, f'{nu}_vtx_wall{vcut}', f"{nu} vertex vs wall; Wall", 5, 1, 6)
    ut.bookHist(h, f'{nu}_vtx_brick{vcut}', f"{nu} vertex vs brick; Brick", 60, 0, 60)
    ut.bookHist(h, f'{nu}_n_prong{vcut}', "N. of charged tracks per vtx; Multiplicity", 50, 0, 50)
    if nu != 'neu':
      ut.bookHist(h, f'{nu}_neu_vtx{vcut}', "N. of neutral secondary vertices; N", 50, 0, 50)
      ut.bookHist(h, f'{nu}_difftheta_TXTY{vcut}', 'Difference between outgoing lepton and incoming nu TXTY;TX;TY', 1000, -0.1, 0.1, 1000, -0.1, 0.1)
      ut.bookHist(h, f'{nu}_lep_Energy{vcut}', f'Energy correlation between outgoing lepton and incoming nu ;E_{nu};E_lep', 240, 0, 4800, 240, 0, 4800)
  for particle in ['mu', 'e', 'hadron', 'hadron_sec']:
    ut.bookHist(h, f'{particle}_energy{vcut}', particle+' Energy;Energy [GeV]', 240, 0, 4800)
    ut.bookHist(h, f'{particle}_TX{vcut}', particle+';TX', 300, -1.5, 1.5)
    ut.bookHist(h, f'{particle}_TY{vcut}', particle+';TY', 300, -1.5, 1.5)
    ut.bookHist(h, f'{particle}_TXTY{vcut}', particle+';TX;TY', 300, -1.5, 1.5, 300, -1.5, 1.5)
    ut.bookHist(h, f'{particle}_PT{vcut}', particle+';PT', 300, 0, 3000)

ntuple = ROOT.TNtuple("cbmsim", "Ntuple of nu",'evID:FLUKA_weight:nu_E:nu_tx:nu_ty:nu_vx:nu_vy:nu_vz:nu_wall:nu_brick:lep_E:lep_tx:lep_ty:n_prong:neu_vtx')


###### EVENT LOOP ##############
for i_event, event in enumerate(sTree):
  if len(sys.argv) == 2:
    if i_event > up_to_ev: break
  elif len(sys.argv) == 3:
    if i_event < from_ev: continue
    if i_event >= to_ev: break
  if i_event%100 == 0: print("Sanity check, current event ", i_event)
  
  is_numu = False
  is_nue = False
  nutrack = event.MCTrack[0]
  leptrack = event.MCTrack[1]
  nupdg = nutrack.GetPdgCode()
  leppdg = leptrack.GetPdgCode()
  if abs(nupdg) == 14 and abs(leppdg) == 13: is_numu = True
  elif abs(nupdg) == 12 and abs(leppdg) == 11: is_nue = True
  else: continue

  #Process identification
  gst.GetEntry(i_event)
  from_charm = gst.FLUKA_weight
  if is_numu:
    nu = 'numu' 
    particle = 'mu'
  elif is_nue:
    nu = 'nue'
    particle = 'e'
  if from_charm ==1: nu = nu+'_c'

  nu_vtx = ROOT.TVector3(nutrack.GetStartX(), nutrack.GetStartY(), nutrack.GetStartZ())
  nu_in_brick, nu_brick_int = getBrickInt(nu_vtx, brickRanges)
  nu_wall_int, nu_brick_int = decodeBrick(nu_brick_int)
  if not nu_in_brick: continue  # excluding neutrinos not interacting in the target
  nu_angle = ROOT.TVector3(nutrack.GetPx()/nutrack.GetPz(), nutrack.GetPy()/nutrack.GetPz(), 1.)
  lep_angle = ROOT.TVector3(leptrack.GetPx()/(leptrack.GetPz()), leptrack.GetPy()/(leptrack.GetPz()), 1.)
  lep_pt = leptrack.GetPt()
  diff_theta = lep_angle-nu_angle
  lep_energy = leptrack.GetEnergy()
  lep_theta = ROOT.TMath.Sqrt(lep_angle.X()**2 + lep_angle.Y()**2)
  for ecut, vcut in vis_cuts.items():
    if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
    h[f'{nu}_energy{vcut}'].Fill(nutrack.GetEnergy())
    h[f'{nu}_TX{vcut}'].Fill(nu_angle.X())
    h[f'{nu}_TY{vcut}'].Fill(nu_angle.Y())
    h[f'{nu}_TXTY{vcut}'].Fill(nu_angle.X(), nu_angle.Y())
    h[f'{nu}_vx{vcut}'].Fill(nutrack.GetStartX())
    h[f'{nu}_vy{vcut}'].Fill(nutrack.GetStartY())
    h[f'{nu}_vz{vcut}'].Fill(nutrack.GetStartZ())
    h[f'{nu}_vxy{vcut}'].Fill(nutrack.GetStartX(), nutrack.GetStartY())
    h[f'{nu}_vxz{vcut}'].Fill(nutrack.GetStartZ(), nutrack.GetStartX())
    h[f'{nu}_vyz{vcut}'].Fill(nutrack.GetStartZ(), nutrack.GetStartY())
    h[f'{nu}_difftheta_TXTY{vcut}'].Fill(diff_theta.X(), diff_theta.Y())
    h[f'{nu}_lep_Energy{vcut}'].Fill(nutrack.GetEnergy(), leptrack.GetEnergy())
    h[f'{nu}_vtx_wall{vcut}'].Fill(nu_wall_int)
    h[f'{nu}_vtx_brick{vcut}'].Fill(nu_brick_int)
    h[f'{particle}_energy{vcut}'].Fill(lep_energy)
    h[f'{particle}_TX{vcut}'].Fill(lep_angle.X())
    h[f'{particle}_TY{vcut}'].Fill(lep_angle.Y())
    h[f'{particle}_TXTY{vcut}'].Fill(lep_angle.X(), lep_angle.Y())
  
  n_prong = {0: 0, 1: 0, 200:0, 500: 0, 1000: 0}
  n_neu_vtx = {0: 0, 1: 0, 200:0, 500: 0, 1000: 0}
  for i_track, track in enumerate(event.MCTrack):
    MotherID = track.GetMotherId()
    Energy = track.GetEnergy()
    PdgCode = track.GetPdgCode()
    if not MyPDG.GetParticle(PdgCode):
        if PdgCode not in failedPDGs: failedPDGs.append(PdgCode)
        continue
    charge = MyPDG.GetParticle(PdgCode).Charge()
    tx = track.GetPx()/track.GetPz()
    ty = track.GetPy()/track.GetPz()
    pt = track.GetPt()
    theta = ROOT.TMath.Sqrt(tx**2 + ty**2)
    if MotherID == 0 and charge != 0:
      for ecut, vcut in vis_cuts.items():
        if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
        h[f'hadron_energy{vcut}'].Fill(Energy)
        h[f'hadron_TX{vcut}'].Fill(tx)
        h[f'hadron_TY{vcut}'].Fill(ty)
        h[f'hadron_TXTY{vcut}'].Fill(tx, ty)
        h[f'hadron_PT{vcut}'].Fill(pt)
        n_prong[ecut] += 1
    if charge == 0:
      n_neu_vtx[ecut] += 1 
      n_prong_sec =  {0: 0, 1: 0, 200:0, 500: 0, 1000: 0}
      for j_track, track2 in enumerate(event.MCTrack):
        MotherID2 = track2.GetMotherId()
        if MotherID2 != i_track: continue   #get daughter from neutral
        neu_vtx = ROOT.TVector3(track2.GetStartX(), track2.GetStartY(), track2.GetStartZ())
        neu_in_brick, neu_brick_int = getBrickInt(neu_vtx, brickRanges)
        neu_wall_int, neu_brick_int = decodeBrick(neu_brick_int)
        if not neu_in_brick: # excluding neutrals not interacting in the target
          n_neu_vtx[ecut] -= 1 
          break 
        PdgCode2 = track2.GetPdgCode()
        if not MyPDG.GetParticle(PdgCode2):
            if PdgCode2 not in failedPDGs: failedPDGs.append(PdgCode2)
            continue
        charge2 = MyPDG.GetParticle(PdgCode2).Charge()
        if charge2 == 0: continue
        Energy2 = track2.GetEnergy()
        tx2 = track2.GetPx()/track2.GetPz()
        ty2 = track2.GetPy()/track2.GetPz()
        pt2 = track2.GetPt()
        theta2 = ROOT.TMath.Sqrt(tx2**2 + ty2**2)
        for ecut, vcut in vis_cuts.items():
          if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
          h[f'hadron_sec_energy{vcut}'].Fill(Energy2)
          h[f'hadron_sec_TX{vcut}'].Fill(tx2)
          h[f'hadron_sec_TY{vcut}'].Fill(ty2)
          h[f'hadron_sec_TXTY{vcut}'].Fill(tx2, ty2)
          h[f'hadron_sec_PT{vcut}'].Fill(pt2)
          n_prong_sec[ecut] +=1
      for ecut, vcut in vis_cuts.items():
        if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
        h[f'neu_energy{vcut}'].Fill(Energy)
        h[f'neu_TX{vcut}'].Fill(tx)
        h[f'neu_TY{vcut}'].Fill(ty)
        h[f'neu_TXTY{vcut}'].Fill(tx, ty)
        h[f'neu_vx{vcut}'].Fill(neu_vtx.X())
        h[f'neu_vy{vcut}'].Fill(neu_vtx.Y())
        h[f'neu_vz{vcut}'].Fill(neu_vtx.Z())
        h[f'neu_vxz{vcut}'].Fill(neu_vtx.Z(), neu_vtx.X())
        h[f'neu_vyz{vcut}'].Fill(neu_vtx.Z(), neu_vtx.Y())
        h[f'neu_vxy{vcut}'].Fill(neu_vtx.X(), neu_vtx.Y())
        h[f'neu_vtx_wall{vcut}'].Fill(neu_wall_int)
        h[f'neu_vtx_brick{vcut}'].Fill(neu_brick_int)
        h[f'neu_n_prong{vcut}'].Fill(n_prong_sec[ecut])
  for ecut, vcut in vis_cuts.items():
    if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
    h[f'{nu}_n_prong{vcut}'].Fill(n_prong[ecut])
    h[f'{nu}_neu_vtx{vcut}'].Fill(n_neu_vtx[ecut])
  ntuple.Fill(i_event, from_charm, nutrack.GetEnergy(), nu_angle.X(), nu_angle.Y(), nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z(), nu_wall_int, nu_brick_int, lep_energy, lep_angle.X(), lep_angle.Y(), n_prong[0], n_neu_vtx[0])
###########################################
print('Arrived at event', i_event-1)
tag = ''
if len(sys.argv) == 2: tag = 'to_evt_'+str(up_to_ev)
elif  len(sys.argv) == 3: tag = 'from_ev_'+str(from_ev)+'_to_ev_'+str(to_ev)
for nu in ['numu', 'nue']:
  outFileName = outPath+'/histo_'+nu+'.'+tag+'.root'
  outFile = ROOT.TFile(outFileName, 'RECREATE')
  for _h in h.values():
    if nu in _h.GetName():
      _h.Write()
  outFile.Write()
  outFile.Close()
outNtupleName = outPath+'/ntuple_.'+tag+'.root'
outNtuple = ROOT.TFile(outNtupleName, 'RECREATE')
ntuple.Write()
outNtuple.Write()
outNtuple.Close() 
print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
print("Generated file ", outFileName)
print("Generated ntuple ", outNtupleName)