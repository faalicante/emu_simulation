import ROOT
import os
import re

def edit_evt_rootrc(file, evID):
  track_line = 'fedra.readCPcut:'
  vertex_line = 'emvertex.vtx.cutvtx:'
  try:
    with open(file, 'r') as f:
      lines = f.readlines()

    with open(file, 'w') as f:
      for line in lines:
        if line.startswith(track_line):
          f.write(f"{track_line} s.eMCEvt=={evID}\n")
        elif line.startswith(vertex_line):
          f.write(f"{vertex_line} s.eMCEvt=={evID}\n")
        else:
          f.write(line)

  except FileNotFoundError:
      print(f"File not found: {file}")
  except Exception as e:
      print(f"An error occurred: {str(e)}")

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
  # if i_event%10 == 0: print("Sanity check, current event ", i_event)
  nutrack = event.MCTrack[0]
  nu_vtx = ROOT.TVector3(nutrack.GetStartX(), nutrack.GetStartY(), nutrack.GetStartZ())
  nu_brick_int= getVolInt(nu_vtx)
  print(f'Event {i_event} brick number {nu_brick_int}')
  if nu_brick_int == None: continue
  brickfolder = 'b'+str(nu_brick_int).zfill(6)
  print(brickfolder)
  # os.system(f'cd {brickfolder}')
  os.chdir(brickfolder)
  os.system('pwd')
  os.system('cp ../track.rootrc ./')
  os.system('cp ../vertex.rootrc ./')
  edit_evt_rootrc('track.rootrc', i_event)
  edit_evt_rootrc('vertex.rootrc', i_event)
  os.system(f'emtra -set={nu_brick_int} -new -v=2')
  os.system(f'emvertex -set={nu_brick_int}.0.0.0 -v=2')
  os.system(f'mv {brickfolder}.0.0.0.trk.root {brickfolder}.0.0.{i_event+1}.trk.root')
  os.system(f'mv {brickfolder}.0.0.0.vtx.root {brickfolder}.0.0.{i_event+1}.vtx.root')
  os.chdir('../')

