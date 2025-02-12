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
parser.add_argument("--nue", dest="nue", required=False, type=int, default=None)
options = parser.parse_args()

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

def scattAngle(nu_angle, lep_angle):
  num = lep_angle.X()*nu_angle.X()+lep_angle.Y()*nu_angle.Y()+1
  denom = ROOT.TMath.Sqrt(lep_angle.X()**2+lep_angle.Y()**2+1)*ROOT.TMath.Sqrt(nu_angle.X()**2+nu_angle.Y()**2+1)
  scatt_angle = ROOT.TMath.ACos(float(num)/float(denom))  ##dot product
  return scatt_angle

def scattPhi(px1, py1, px2, py2):
  phi1 = ROOT.TMath.ATan2(py1, px1)
  phi2 = ROOT.TMath.ATan2(py2, px2)
  scatt_phi = ROOT.TMath.Abs(phi1-phi2)
  if scatt_phi > ROOT.TMath.Pi(): scatt_phi = 2*ROOT.TMath.Pi() - scatt_phi
  return scatt_phi

def getVolInt(nu_vtx):
  nodeInt = ROOT.gGeoManager.FindNode(nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z())
  pathInt = ROOT.gGeoManager.GetPath()

  wall_path = re.search(r"/Wall_(\d+)", pathInt)
  row_path = re.search(r"/Row_(\d+)", pathInt)
  brick_path = re.search(r"/Brick_(\d+)", pathInt)

  wall_number = int(wall_path.group(1)) if wall_path else None
  row_number = int(row_path.group(1)) if row_path else None
  brick_number = int(brick_path.group(1)) if brick_path else None
  if brick_number != None: wall, brick = decodeVolInt(wall_number, row_number, brick_number)
  else: return None, None, pathInt
  return wall, brick, pathInt

def decodeVolInt(wall_number, row_number, brick_number):
  wall = wall_number + 1
  column = brick_number
  row = row_number
  brick_map = {(0, 0): 2, (0, 1): 1, (1, 0): 4, (1, 1): 3}
  brick = brick_map[(row, column)]
  brick = wall*10 + brick
  return wall, brick


numu = options.numu
nue = options.nue
from_evt = options.from_evt
to_evt = options.to_evt

outPath = '/afs/cern.ch/work/f/falicant/public/emu_simulation/out2'

if numu:
  pathSim = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc_eff10_smear0'
  nuf = 'numu'
  lep = 'mu'
  flag = 1
elif nue:
  pathSim = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc_eff10_smear0'
  nuf = 'nue'
  lep = 'e'
  flag = 2
else:
  print('No input')
  sys.exit(1)

geoFile =  pathSim + '/geofile_full.Genie-TGeant4.root'
# geo = SndlhcGeo.GeoInterface(geoFile)
ROOT.TGeoManager.Import(geoFile)

simName = '/inECC_sndLHC.Genie-TGeant4.root'
simFile = ROOT.TFile.Open(pathSim+simName)
sTree = simFile.cbmsim
# gst = simFile.gst
nentries = sTree.GetEntries()
print('There are ', nentries, 'entries')

wallRanges = getWallRanges()
brickRanges = getBrickRanges()

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()


tag = 'from_evt_'+str(from_evt)+'_to_evt_'+str(to_evt)
outNtupleName = outPath+f'/ntuple_{nuf}_b11.{tag}.root'
outNtuple = ROOT.TFile.Open(outNtupleName, 'RECREATE')
ntuple = ROOT.TNtuple("cbmsim", "Ntuple of nu",'evID:flag:nu_E:nu_vx:nu_vy:nu_vz:nu_wall:nu_brick:lep_E:lep_tx:lep_ty:n_prong:neu_vtx')
ntuple.SetDirectory(outNtuple)

h               = {}
e_cuts          = {200:'_e200', 500:'_e500', 1000:'_e1000'}
angle_cuts      = {0:'', 1:'_angle1'}
vis_cuts        = {**angle_cuts, **e_cuts}

for nu in [nuf, 'neu']:
  ut.bookHist(h, f'{nu}_vx_all', f";{nu} vtx X[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
  ut.bookHist(h, f'{nu}_vy_all', f";{nu} vtx Y[cm]", int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
  ut.bookHist(h, f'{nu}_vz_all', f";{nu} vtx Z[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20)
  ut.bookHist(h, f'{nu}_lost', f";{nu} vtx outside of target", 2,0,1)
  ut.bookHist(h, f'{nu}_energy', f';E_{nu}', 200, 0, 4000)
  ut.bookHist(h, f'{nu}_TX', ';TX', 2000, -1, 1)
  ut.bookHist(h, f'{nu}_TY', ';TY', 2000, -1, 1)
  ut.bookHist(h, f'{nu}_TXTY', ';TX;TY', 2000, -1, 1, 2000, -1, 1)
  ut.bookHist(h, f'{nu}_vx', f";{nu} vtx X[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
  ut.bookHist(h, f'{nu}_vy', f";{nu} vtx Y[cm]", int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
  ut.bookHist(h, f'{nu}_vz', f";{nu} vtx Z[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20)
  ut.bookHist(h, f'{nu}_vxy', f";{nu} vtx X[cm];{nu} vtx Y[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
  ut.bookHist(h, f'{nu}_vxz', f";{nu} vtx Z[cm];{nu} vtx X[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
  ut.bookHist(h, f'{nu}_vyz', f";{nu} vtx Z[cm];{nu} vtx Y[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
  ut.bookHist(h, f'{nu}_vtx_wall', f"{nu} vertex vs wall; Wall", 5, 1, 6)
  ut.bookHist(h, f'{nu}_vtx_brick', f"{nu} vertex vs brick; Brick", 60, 0, 60)
  for vcut in vis_cuts.values():
    ut.bookHist(h, f'{nu}_n_prong{vcut}', "N. of charged tracks per vtx; Multiplicity", 50, 0, 50)
ut.bookHist(h, f'{nuf}_neu_vtx', "N. of neutral secondary vertices; N", 100, 0, 100)
ut.bookHist(h, f'{nuf}_scatt_ang', 'Scattering angle between outgoing lepton and incoming nu;#theta;#phi', 320, 3.2, 3.2, 320, 3.2, 3.2)
ut.bookHist(h, f'{nuf}_lep_Energy', f'Energy correlation between outgoing lepton and incoming nu ;E_{nu};E_lep', 200, 0, 4000, 200, 0, 4000)
for particle in [lep, 'hadron', 'hadron_sec']:
  for vcut in vis_cuts.values():
    ut.bookHist(h, f'{particle}_energy{vcut}', particle+' Energy;Energy [GeV]', 200, 0, 4000)
    ut.bookHist(h, f'{particle}_PT{vcut}', particle+';PT', 1000, 0, 100)
  ut.bookHist(h, f'{particle}_TX', particle+';TX', 2000, -1, 1)
  ut.bookHist(h, f'{particle}_TY', particle+';TY', 2000, -1, 1)
  ut.bookHist(h, f'{particle}_TXTY', particle+';TX;TY', 2000, -1, 1, 2000, -1, 1)


###### EVENT LOOP ##############
for i_event, event in enumerate(sTree):
  if i_event < from_evt: continue
  if i_event >= to_evt: break
  if i_event%100 == 0: print("Sanity check, current event ", i_event)

  nutrack = event.MCTrack[0]
  leptrack = event.MCTrack[1]
  nupdg = nutrack.GetPdgCode()
  leppdg = leptrack.GetPdgCode()

  nu_vtx = ROOT.TVector3(nutrack.GetStartX(), nutrack.GetStartY(), nutrack.GetStartZ())
  nu_wall_int, nu_brick_int, vol_path_int = getVolInt(nu_vtx)
  h[f'{nu}_vx_all'].Fill(nu_vtx.X())
  h[f'{nu}_vy_all'].Fill(nu_vtx.Y())
  h[f'{nu}_vz_all'].Fill(nu_vtx.Z())
  # nu_in_brick, nu_brick_int = getBrickInt(nu_vtx, brickRanges)
  # nu_wall_int, nu_brick_int = decodeBrick(nu_brick_int)
  if nu_brick_int == None: # excluding neutrinos not interacting in the target
    h[f'{nu}_lost'].Fill(0)
    print(f"Vtx {i_event} not in brick", nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z())
    print("Vtx in volume", vol_path_int)
    continue  
  if nu_brick_int != 11: continue
  nu_angle = ROOT.TVector3(nutrack.GetPx()/nutrack.GetPz(), nutrack.GetPy()/nutrack.GetPz(), 1.)
  lep_angle = ROOT.TVector3(leptrack.GetPx()/(leptrack.GetPz()), leptrack.GetPy()/(leptrack.GetPz()), 1.)
  lep_pt = leptrack.GetPt()
  scatt_angle = scattAngle(nu_angle, lep_angle)
  scatt_phi = scattPhi(nutrack.GetPx(), nutrack.GetPy(), leptrack.GetPx(), leptrack.GetPy())
  nu_energy = nutrack.GetEnergy()
  lep_energy = leptrack.GetEnergy()
  lep_theta = ROOT.TMath.Sqrt(lep_angle.X()**2 + lep_angle.Y()**2)
  h[f'{nuf}_energy'].Fill(nu_energy)
  h[f'{nuf}_TX'].Fill(nu_angle.X())
  h[f'{nuf}_TY'].Fill(nu_angle.Y())
  h[f'{nuf}_TXTY'].Fill(nu_angle.X(), nu_angle.Y())
  h[f'{nuf}_vx'].Fill(nu_vtx.X())
  h[f'{nuf}_vy'].Fill(nu_vtx.Y())
  h[f'{nuf}_vz'].Fill(nu_vtx.Z())
  h[f'{nuf}_vxy'].Fill(nu_vtx.X(), nu_vtx.Y())
  h[f'{nuf}_vxz'].Fill(nu_vtx.Z(), nu_vtx.X())
  h[f'{nuf}_vyz'].Fill(nu_vtx.Z(), nu_vtx.Y())
  h[f'{nuf}_scatt_ang'].Fill(scatt_angle, scatt_phi)
  h[f'{nuf}_lep_Energy'].Fill(nu_energy, lep_energy)
  h[f'{nuf}_vtx_wall'].Fill(nu_wall_int)
  h[f'{nuf}_vtx_brick'].Fill(nu_brick_int)
  h[f'{lep}_TX'].Fill(lep_angle.X())
  h[f'{lep}_TY'].Fill(lep_angle.Y())
  h[f'{lep}_TXTY'].Fill(lep_angle.X(), lep_angle.Y())
  for ecut, vcut in vis_cuts.items():
    if ecut and (lep_theta>1 or lep_energy < ecut/1e3): continue
    h[f'{lep}_energy{vcut}'].Fill(lep_energy)
    h[f'{lep}_PT{vcut}'].Fill(lep_pt)
  
  n_prong = {0:0, 1:0, 200:0, 500:0, 1000:0}
  n_neu_vtx = 0
  for i_track, track in enumerate(event.MCTrack):
    if i_track < 1 : continue #skip neutrino
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
      h[f'hadron_TX'].Fill(tx)
      h[f'hadron_TY'].Fill(ty)
      h[f'hadron_TXTY'].Fill(tx, ty)
      for ecut, vcut in vis_cuts.items():
        if ecut and (theta>1 or Energy < ecut/1e3): continue
        n_prong[ecut] += 1
        if i_track == 1: continue #skip lepton
        h[f'hadron_energy{vcut}'].Fill(Energy)
        h[f'hadron_PT{vcut}'].Fill(pt)
    if charge == 0:
      n_prong_sec = {0:0, 1:0, 200:0, 500:0, 1000:0}
      for j_track, track2 in enumerate(event.MCTrack):
        MotherID2 = track2.GetMotherId()
        if MotherID2 != i_track: continue   #get daughter from neutral
        neu_vtx = ROOT.TVector3(track2.GetStartX(), track2.GetStartY(), track2.GetStartZ())
        neu_in_brick, neu_brick_int = getBrickInt(neu_vtx, brickRanges)
        neu_wall_int, neu_brick_int = decodeBrick(neu_brick_int)
        if not neu_in_brick: break # excluding neutrals not interacting in the target
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
        h[f'hadron_sec_TX'].Fill(tx2)
        h[f'hadron_sec_TY'].Fill(ty2)
        h[f'hadron_sec_TXTY'].Fill(tx2, ty2)
        for ecut, vcut in vis_cuts.items():
          if ecut and (theta2>1 or Energy2 < ecut/1e3): continue
          h[f'hadron_sec_energy{vcut}'].Fill(Energy2)
          h[f'hadron_sec_PT{vcut}'].Fill(pt2)
          n_prong_sec[ecut] +=1
      for ecut, vcut in vis_cuts.items():
        if not n_prong_sec[ecut]: continue
        h[f'neu_n_prong{vcut}'].Fill(n_prong_sec[ecut])
      if not n_prong_sec[1]: continue
      n_neu_vtx += 1
      h[f'neu_energy'].Fill(Energy)
      h[f'neu_TX'].Fill(tx)
      h[f'neu_TY'].Fill(ty)
      h[f'neu_TXTY'].Fill(tx, ty)
      h[f'neu_vx'].Fill(neu_vtx.X())
      h[f'neu_vy'].Fill(neu_vtx.Y())
      h[f'neu_vz'].Fill(neu_vtx.Z())
      h[f'neu_vxz'].Fill(neu_vtx.Z(), neu_vtx.X())
      h[f'neu_vyz'].Fill(neu_vtx.Z(), neu_vtx.Y())
      h[f'neu_vxy'].Fill(neu_vtx.X(), neu_vtx.Y())
      h[f'neu_vtx_wall'].Fill(neu_wall_int)
      h[f'neu_vtx_brick'].Fill(neu_brick_int)
  
  for ecut, vcut in vis_cuts.items():
    h[f'{nuf}_n_prong{vcut}'].Fill(n_prong[ecut])
  h[f'{nuf}_neu_vtx'].Fill(n_neu_vtx)
  ntuple.Fill(i_event, flag, nu_energy, nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z(), nu_wall_int, nu_brick_int, lep_energy, lep_angle.X(), lep_angle.Y(), n_prong[0], n_neu_vtx)
###########################################
print('Arrived at event', i_event-1)

outNtuple.cd()
ntuple.Write()
outNtuple.Write()
outNtuple.Close() 

outFileName = outPath+f'/histo_{nuf}.{tag}.root'
outFile = ROOT.TFile(outFileName, 'RECREATE')
for _h in h.values():
  _h.Write()
outFile.Write()
outFile.Close()

print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
print("Generated file ", outFileName)
print("Generated ntuple ", outNtupleName)