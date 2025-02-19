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
parser.add_argument("--K", dest="K", required=False, type=bool, default=None)
parser.add_argument("--neu", dest="neu", required=False, type=bool, default=None)
parser.add_argument("--part", dest="part", required=True, type=int, default=None)
parser.add_argument("--enbin", dest="enbin", required=True, type=int, default=None)
options = parser.parse_args()

K = options.K
neu = options.neu
part = options.part
enbin = options.enbin

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

outPath = '/afs/cern.ch/work/f/falicant/public/emu_simulation/out_neutrals'

prepath = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/neutrals/'
energyBins = ['_10_50', '_50_100', '_100_200']

if K:
  pathSim = prepath+'K'+energyBins[enbin]+'/Ntuples/Ntuples/'+str(part)
  simName = '/sndLHC.PG_130-TGeant4.root'
  geoFile =  pathSim + '/geofile_full.PG_130-TGeant4.root'
  flavour = 'K'
  flag = 10+enbin
elif neu:
  pathSim = prepath+'neu'+energyBins[enbin]+'/Ntuples/Ntuples/'+str(part)
  simName = '/sndLHC.PG_2112-TGeant4.root'
  geoFile =  pathSim + '/geofile_full.PG_2112-TGeant4.root'
  flavour = 'neu'
  flag = 20+enbin
else:
  print('No input')
  sys.exit(1)

# geo = SndlhcGeo.GeoInterface(geoFile)
ROOT.TGeoManager.Import(geoFile)
simFile = ROOT.TFile.Open(pathSim+simName)
sTree = simFile.cbmsim
nentries = sTree.GetEntries()
print('There are ', nentries, 'entries')

wallRanges = getWallRanges()
# brickRanges = getBrickRanges()

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()

tag = flavour+energyBins[enbin]+'_'+str(part)
outNtupleName = outPath+f'/ntuple.{tag}.root'
outNtuple = ROOT.TFile.Open(outNtupleName, 'RECREATE')
ntuple = ROOT.TNtuple("cbmsim", "Ntuple of nu",'evID:flag:energy:vx:vy:vz:wall:brick:n_prong')
ntuple.SetDirectory(outNtuple)

h               = {}
e_cuts          = {200:'_e200', 500:'_e500', 1000:'_e1000'}
angle_cuts      = {0:'', 1:'_angle1'}
vis_cuts        = {**angle_cuts, **e_cuts}

particle = 'hadron'

ut.bookHist(h, f'{flavour}_vx_all', f";{flavour} vtx X[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
ut.bookHist(h, f'{flavour}_vy_all', f";{flavour} vtx Y[cm]", int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
ut.bookHist(h, f'{flavour}_vz_all', f";{flavour} vtx Z[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20)
ut.bookHist(h, f'{flavour}_lost', f";{flavour} vtx outside of target", 2,0,1)
ut.bookHist(h, f'{flavour}_energy', f';E_{flavour}', 200, 0, 4000)
ut.bookHist(h, f'{flavour}_TX', ';TX', 2000, -1, 1)
ut.bookHist(h, f'{flavour}_TY', ';TY', 2000, -1, 1)
ut.bookHist(h, f'{flavour}_TXTY', ';TX;TY', 2000, -1, 1, 2000, -1, 1)
ut.bookHist(h, f'{flavour}_vx', f";{flavour} vtx X[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
ut.bookHist(h, f'{flavour}_vy', f";{flavour} vtx Y[cm]", int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
ut.bookHist(h, f'{flavour}_vz', f";{flavour} vtx Z[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20)
ut.bookHist(h, f'{flavour}_vxy', f";{flavour} vtx X[cm];{flavour} vtx Y[cm]", int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
ut.bookHist(h, f'{flavour}_vxz', f";{flavour} vtx Z[cm];{flavour} vtx X[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['X'][4][1]-wallRanges['X'][0][0]+40)*10, wallRanges['X'][0][0]-20, wallRanges['X'][4][1]+20)
ut.bookHist(h, f'{flavour}_vyz', f";{flavour} vtx Z[cm];{flavour} vtx Y[cm]", int(wallRanges['Z'][4][1]-wallRanges['Z'][0][0]+40)*10, wallRanges['Z'][0][0]-20, wallRanges['Z'][4][1]+20, int(wallRanges['Y'][4][1]-wallRanges['Y'][0][0]+40)*10, wallRanges['Y'][0][0]-20, wallRanges['Y'][4][1]+20)
ut.bookHist(h, f'{flavour}_vtx_wall', f"{flavour} vertex vs wall; Wall", 5, 1, 6)
ut.bookHist(h, f'{flavour}_vtx_brick', f"{flavour} vertex vs brick; Brick", 60, 0, 60)
for vcut in vis_cuts.values():
  ut.bookHist(h, f'{flavour}_n_prong{vcut}', "N. of charged tracks per vtx; Multiplicity", 50, 0, 50)
  ut.bookHist(h, f'{particle}_energy{vcut}', particle+' Energy;Energy [GeV]', 200, 0, 4000)
  ut.bookHist(h, f'{particle}_PT{vcut}', particle+';PT', 1000, 0, 100)
  ut.bookHist(h, f'{particle}_TX{vcut}', particle+';TX', 2000, -1, 1)
  ut.bookHist(h, f'{particle}_TY{vcut}', particle+';TY', 2000, -1, 1)
  ut.bookHist(h, f'{particle}_TXTY{vcut}', particle+';TX;TY', 2000, -1, 1, 2000, -1, 1)

###### EVENT LOOP ##############
for i_event, event in enumerate(sTree):
  if i_event%10 == 0: print("Sanity check, current event ", i_event)

  neutrack = event.MCTrack[0]
  inttrack = event.MCTrack[1]
  neu_vtx = ROOT.TVector3(inttrack.GetStartX(), inttrack.GetStartY(), inttrack.GetStartZ())
  neu_wall_int, neu_brick_int, vol_path_int = getVolInt(neu_vtx)
  h[f'{flavour}_vx_all'].Fill(neu_vtx.X())
  h[f'{flavour}_vy_all'].Fill(neu_vtx.Y())
  h[f'{flavour}_vz_all'].Fill(neu_vtx.Z())
  if neu_brick_int == None: # excluding neutrinos not interacting in the target
    h[f'{flavour}_lost'].Fill(0)
    print(f"Vtx {i_event} not in brick", vol_path_int)
    continue  
  neu_angle = ROOT.TVector3(neutrack.GetPx()/neutrack.GetPz(), neutrack.GetPy()/neutrack.GetPz(), 1.)
  neu_energy = neutrack.GetEnergy()
  h[f'{flavour}_energy'].Fill(neu_energy)
  h[f'{flavour}_TX'].Fill(neu_angle.X())
  h[f'{flavour}_TY'].Fill(neu_angle.Y())
  h[f'{flavour}_TXTY'].Fill(neu_angle.X(), neu_angle.Y())
  h[f'{flavour}_vx'].Fill(neu_vtx.X())
  h[f'{flavour}_vy'].Fill(neu_vtx.Y())
  h[f'{flavour}_vz'].Fill(neu_vtx.Z())
  h[f'{flavour}_vxy'].Fill(neu_vtx.X(), neu_vtx.Y())
  h[f'{flavour}_vxz'].Fill(neu_vtx.Z(), neu_vtx.X())
  h[f'{flavour}_vyz'].Fill(neu_vtx.Z(), neu_vtx.Y())
  h[f'{flavour}_vtx_wall'].Fill(neu_wall_int)
  h[f'{flavour}_vtx_brick'].Fill(neu_brick_int)
  
  n_prong = {0:0, 1:0, 200:0, 500:0, 1000:0}
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
        if ecut and (theta>1 or Energy < ecut/1e3): continue
        n_prong[ecut] += 1
        h[f'hadron_energy{vcut}'].Fill(Energy)
        h[f'hadron_TX{vcut}'].Fill(tx)
        h[f'hadron_TY{vcut}'].Fill(ty)
        h[f'hadron_TXTY{vcut}'].Fill(tx, ty)
        h[f'hadron_PT{vcut}'].Fill(pt)
  for ecut, vcut in vis_cuts.items():
    h[f'{flavour}_n_prong{vcut}'].Fill(n_prong[ecut])
  ntuple.Fill(i_event, flag, neu_energy, neu_vtx.X(), neu_vtx.Y(), neu_vtx.Z(), neu_wall_int, neu_brick_int, n_prong[0])
###########################################
print('Arrived at event', i_event)

outNtuple.cd()
ntuple.Write()
outNtuple.Write()
outNtuple.Close() 

outFileName = outPath+f'/histo.{tag}.root'
outFile = ROOT.TFile(outFileName, 'RECREATE')
for _h in h.values():
  _h.Write()
outFile.Write()
outFile.Close()

print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
print("Generated file ", outFileName)
print("Generated ntuple ", outNtupleName)