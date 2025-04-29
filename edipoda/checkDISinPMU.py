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
    DISevts = {}
    for ievent, event in enumerate(cbmsim):
        if ievent%1000==0: print("event n. ", ievent)
        for itrack, track in enumerate(event.MCTrack):
            if track.GetMotherId()==0 and track.GetProcID() in DIS_procIDs:
                DISevts.setdefault(ievent, []).append(itrack)
    return DISevts


brickid = options.brick
simPath = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1E5/muon1E5_eff9_smear1'
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

cbmsim_file = ROOT.TFile.Open(simPath+'/sndLHC.Ntuple-TGeant4-1E5cm2.root')
cbmsim = cbmsim_file.cbmsim
DISevts = prepareList(cbmsim)
print(DISevts)

vtx_arr_dis = ROOT.TObjArray()
for ivtx, vtx in enumerate(vertices):
    print("Vertex n.", ivtx)
    foundDIS = False
    vID = vtx.ID()
    ntrks = vtx.N()
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        nseg = track.N()
        trackPDG, trackEvt, trackID, trackMother = getMajorityFeatures(track)
        if trackEvt in DISevts.keys() and trackID in DISevts[trackEvt]:
            foundDIS = True
            print(trackEvt, trackID, vID)
            print(DISevts[trackEvt])
            break
    if foundDIS: vtx_arr_dis.Add(vtx)

dproc.MakeVertexTree(vtx_arr_dis, simPath+f"/b0000{brickid}/DIS_in_PMU.vtx.root")
vtx_arr_dis.Delete()  # Explicitly delete to avoid memory leaks
