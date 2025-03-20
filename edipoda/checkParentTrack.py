import ROOT
import fedrarootlogon
from array import array
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b", dest="brick", required=True, type=int)
options = parser.parse_args()

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

def getMajorityFeatures(track):
    nseg = track.N()
    DictSegPDG = {}
    DictSegEvt = {}
    DictSegID = {}
    DictSegMother = {}
    for iseg in range(nseg):
        seg = track.GetSegment(iseg)
        sPDG = seg.Vid(0)
        sEvt = seg.MCEvt()
        sID = seg.MCTrack()
        sMother = seg.Aid(0)
        AddToDict(DictSegPDG, sPDG)
        AddToDict(DictSegEvt, sEvt)
        AddToDict(DictSegID, sID)
        AddToDict(DictSegMother, sMother)
    trackPDG = max(DictSegPDG, key=DictSegPDG.get)
    trackEvt = max(DictSegEvt, key=DictSegEvt.get)
    trackID = max(DictSegID, key=DictSegID.get)
    trackMother = max(DictSegMother, key=DictSegMother.get)
    return trackPDG, trackEvt, trackID, trackMother

def prepareList(cbmsim):
    DIS_procIDs = [23,13,25,24,26,27,46]
    nu_list = []
    for ievt, event in enumerate(cbmsim):
        if ievt%100==0: print(f"Event {ievt}")
        nu_vtx = ROOT.TVector3(event.MCTrack[0].GetStartX(), event.MCTrack[0].GetStartY(), event.MCTrack[0].GetStartZ())
        nu_list.append(nu_vtx)
    return nu_list

def findVertex(vtx, dis_list):
    max_xy=0.0003 #3um
    max_z=0.1350 #1350um
    vtx_conv = convertVertex(vtx, brickID)
    vtx_g = ROOT.TVector3(vtx_conv[0], vtx_conv[1], vtx_conv[2])
    closestEvent = None
    for ievt, nu_vtx in enumerate(dis_list):
        dist_xy = ROOT.TMath.Sqrt((nu_vtx.X()-vtx_g.X())**2 + (nu_vtx.Y()-vtx_g.Y())**2)
        dist_z = ROOT.TMath.Abs(nu_vtx.Z()-vtx_g.Z())
        if dist_z < max_z:
            if dist_xy < max_xy:
                max_xy = dist_xy
                closestEvent = ievt
    return closestEvent

brickid = options.brick
simPath = '/eos/experiment/sndlhc/users/falicant/simulations/muonDIS'
vtx_file = simPath+f'/b0000{brickid}/b0000{brickid}.0.0.0.vtx.root'
print("opening file: ",vtx_file)
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
scancond.SetSigma0(1.5,1.5,0.0015,0.0015) #change sigma0
scancond.SetDegrad(3) #change angular degradation
gAli.SetScanCond(scancond)
vertexrec = ROOT.EdbVertexRec()
vertexrec.SetPVRec(gAli)
vertexrec.eDZmax=3000.
vertexrec.eProbMin=0.001
vertexrec.eImpMax=3.5
vertexrec.eUseMom=False
vertexrec.eUseSegPar=False
vertexrec.eQualityMode=0
proc = ROOT.EdbDataProc()
dproc.ReadVertexTree(vertexrec, vtx_file, "flag==0&&n>2")
vertices = gAli.eVTX

trk_file = simPath+f'/b0000{brickid}/b0000{brickid}.0.0.0.trk.root'
print("opening file: ",trk_file)
dproc2 = ROOT.EdbDataProc()
gAli2 = dproc2.PVR()
scancond2 = ROOT.EdbScanCond()
scancond2.SetSigma0(1.5,1.5,0.0015,0.0015) #change sigma0
scancond2.SetDegrad(3) #change angular degradation
gAli2.SetScanCond(scancond2)
dproc2.ReadTracksTree(gAli2, trk_file, "nseg>2")
tracks = gAli2.eTracks

ftrk_file_name = simPath+f'/b0000{brickid}/found_tracks.root'
ftrk_file = ROOT.TFile.Open(ftrk_file_name)
parents = ftrk_file.tracks

DIS_vIDs = list()
DISevts = list()
foundDIS = 0
for ivtx, vtx in enumerate(vertices):
    # print("Vertex n.", ivtx)
    trackFound = False
    vID = vtx.ID()
    ntrks = vtx.N()
    DictVtxEvt = {}
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        nseg = track.N()
        trackPDG, trackEvt, trackID, trackMother = getMajorityFeatures(track)
        AddToDict(DictVtxEvt, trackEvt)
    vtxEvt = max(DictVtxEvt, key=DictVtxEvt.get)
    for ftrack in parents:
        if not ftrack.chosen: continue
        if ftrack.vid == vID: 
            parentID = int(ftrack.tid)
            trackFound = True
            # print(f"Vtx {vID}, parent {parentID}")
            break
    if not trackFound: continue
    for trackp in tracks:
        if trackp.ID() != parentID: continue
        trackPDGp, trackEvtp, trackIDp, trackMotherp = getMajorityFeatures(trackp)
        break
    if trackEvtp == vtxEvt:
        foundDIS += 1
    else: 
        print("FALSE PARENT!!!  TrackID:", parentID)
        parents.Scan("vid:tid:nseg:firstp:r2:dz",f"vid=={vID}&&tid=={parentID}&&chosen==1")

print('Total vtx:', ivtx+1)
print('Total tracks found:', parents.GetEntries("chosen==1"))
print('Total true parents:', foundDIS)
