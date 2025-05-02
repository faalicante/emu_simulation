import ROOT
import os
import re

def add_to_dat(event, cell, xpos, ypos):
  with open('nu_evt_cells.dat', 'a') as f:
    f.write(f'{event}, {cell}, {xpos}, {ypos}\n')

def decodeVolInt(wall_number, row_number, brick_number):
  wall = wall_number + 1
  column = brick_number
  row = row_number
  brick_map = {(0, 0): 2, (0, 1): 1, (1, 0): 4, (1, 1): 3}
  brick = brick_map[(row, column)]
  brick = wall*10 + brick
  return brick

def getVolInt(nu_vtx):
  nodeInt = ROOT.gGeoManager.FindNode(nu_vtx.X(), nu_vtx.Y(), nu_vtx.Z())
  pathInt = ROOT.gGeoManager.GetPath()

  wall_path = re.search(r"/Wall_(\d+)", pathInt)
  row_path = re.search(r"/Row_(\d+)", pathInt)
  brick_path = re.search(r"/Brick_(\d+)", pathInt)

  wall_number = int(wall_path.group(1)) if wall_path else None
  row_number = int(row_path.group(1)) if row_path else None
  brick_number = int(brick_path.group(1)) if brick_path else None
  if brick_number != None: brick = decodeVolInt(wall_number, row_number, brick_number)
  else: return None
  return brick

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--from_evt", dest="from_evt", required=False, type=int, default=0)
parser.add_argument("--to_evt", dest="to_evt", required=True, type=int, default=None)
options = parser.parse_args()

from_evt = options.from_evt
to_evt = options.to_evt

geoFile = 'geofile_full.Genie-TGeant4.root'
ROOT.TGeoManager.Import(geoFile)

simName = 'inECC_sndLHC.Genie-TGeant4.root'
simFile = ROOT.TFile.Open(simName)
events = simFile.cbmsim

for i_event, event in enumerate(events):
  if i_event < from_evt: continue
  if i_event >= to_evt: break
  print(f'Processing event {i_event}')
  nutrack = event.MCTrack[0]
  nu_vtx = ROOT.TVector3(nutrack.GetStartX(), nutrack.GetStartY(), nutrack.GetStartZ())
  nu_brick_int= getVolInt(nu_vtx)
  if nu_brick_int != 21: continue
  xem = nu_vtx.X() + 27.3
  yem = nu_vtx.Y() - 15.8
  xcell = round(xem-0.5)
  ycell = round(yem)+1
  mucell = ycell*21+xcell
  xpos = round((xem + 20) * 1E4,0)
  ypos = round(yem * 1E4,0)
  add_to_dat(i_event, mucell, xpos, ypos)
