import ROOT
import numpy as np
import fedrarootlogon
from time import time
import os.path
from array import array


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-p", dest="part", required=True, default=None, type=int)
parser.add_argument("-b", dest="brickID", required=True, default=None, type=int)
options = parser.parse_args()

brickID = options.brickID
npart = options.part
from_plate = 60
zmin = -77585.00
zmax = 0.0

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

parts = ['K_10_50', 'K_50_100', 'K_100_200', 'neu_10_50', 'neu_50_100', 'neu_100_200']

##WEIGHT USED FOR DATA IN RUN1 BRICK21
events_b21 = [181, 206, 203, 195, 189, 214]
expected_b21 = [42.840602, 5.7650435, 2.4583636, 24.623528, 2.1840258, 0.40455220]
events_b31 = [118, 93, 105, 110, 100, 99]
expected_b31 = [27.929232, 2.6026653, 1.2715674, 13.890195, 1.1555692, 0.18715265]
w = expected_b21[npart]/events_b21[npart]

path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/neutrals/in_muon1.3E5/'
out_dir = path
out_name = f'/vertex_neutral_{npart}_{brickID}.root'

N=40
start_time = time()

dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
scancond.SetSigma0(0.5,0.5,0.0025,0.0025)
scancond.SetDegrad(5)
gAli.SetScanCond(scancond)
vertexrec = ROOT.EdbVertexRec()
vertexrec.SetPVRec(gAli)
vertexrec.eDZmax=3000.
vertexrec.eProbMin=0.001
vertexrec.eImpMax=3.5
vertexrec.eUseMom=False
vertexrec.eUseSegPar=False
vertexrec.eQualityMode=0

vertextot = ROOT.EdbPVRec()
vtx_cut = f'n>2&&(flag==0||flag==3)&&vz>{zmin}&&vz<{zmax}'

rec_vtx = 0
neu_vtx = 0
for i in range(1000):
    vtx_file = path+parts[npart]+f'/b{brickID:06d}/b{brickID:06d}.0.0.{i+1}.vtx.root'
    print("opening file: ", vtx_file)
    if not os.path.isfile(vtx_file):
        print("Vertex file not found, ", vtx_file)
        continue
    dproc.ReadVertexTree(vertexrec, vtx_file, vtx_cut)
    rec_vtx +=1

for vtx in gAli.eVTX:
    DictFlag = {}
    DictMoth = {}
    for itrk in range(vtx.N()):
        track = vtx.GetTrack(itrk)
        nseg = track.N()
        for iseg in range(nseg):
            seg = track.GetSegment(iseg)
            signal = seg.Flag()
            mother = seg.Aid(0)
            AddToDict(DictFlag, signal)
            AddToDict(DictMoth, mother)
    if max(DictFlag, key=DictFlag.get) == 1 and max(DictMoth, key=DictMoth.get) == 0:
        vertextot.AddVertex(vtx)
        neu_vtx+=1

#save vertices in root file
outputFile = ROOT.TFile(out_dir+out_name,"RECREATE")	
outputTree = ROOT.TTree("vtx","Tree of vertices")
outputTree.SetDirectory(outputFile)

_brickID = array('i', [0])
_npart = array('i', [0])
_vID = array('i', [0])
_flag = array('i', [0])
_vx = array('f', [0])
_vy = array('f', [0])
_vz = array('f', [0])
_ntrks = array('i', [0])
_nseg = array('i', N*[0])
_npl = array('i', N*[0])
_fillfact_t = array('f', N*[0])
_fillfact = array('f', [0])
_meanIP = array('f', [0])
_ip = array('f', N*[0])
_prob = array('f', [0])
_maxaperture = array('f', [0])
_meanaperture = array('f', [0])
_dphi = array('f', N*[0])
_maxdphi = array('f', [0])
_magdphi = array('f', [0])
_meanphi = array('f', [0])
_plate = array('i', [0])
_nfirst = array('i', N*[0])
_nlast = array('i', N*[0])
_tTX = array('f', N*[0])
_tTY = array('f', N*[0])
_tLastX = array('f', N*[0])
_tLastY = array('f', N*[0])
_tLastZ = array('f', N*[0])
_tID = array('i', N*[0])
_tPDG = array('i', N*[0])
_tEvt = array('i', N*[0])
_tMot = array('i', N*[0])
_trPDG = array('f', [0])
_motherdPhi = array('f', [0])
_signal = array('i', [0])
_weight = array('f', [0])
_clEvt = array('i', [0])

outputTree.Branch("brickID", _brickID, "brickID/I")
outputTree.Branch("naprt", _npart, "npart/I")
outputTree.Branch("vID", _vID, "vID/I")
outputTree.Branch("flag", _flag, "flag/I")
outputTree.Branch("vx", _vx, "vx/F")
outputTree.Branch("vy", _vy, "vy/F")
outputTree.Branch("vz", _vz, "vz/F")
outputTree.Branch("ntrks", _ntrks, "ntrks/I")
outputTree.Branch("nseg", _nseg, "nseg[ntrks]/I")
outputTree.Branch("npl", _npl, "npl[ntrks]/I")
outputTree.Branch("fillfact", _fillfact, "fillfact/F")
outputTree.Branch("fillfact_t", _fillfact_t, "fillfact_t[ntrks]/F")
outputTree.Branch("meanIP", _meanIP, "meanIP/F")
outputTree.Branch("ip", _ip, "ip[ntrks]/F")
outputTree.Branch("tTX", _tTX, "tTX[ntrks]/F")
outputTree.Branch("tTY", _tTY, "tTY[ntrks]/F")
outputTree.Branch("tID", _tID, "tID[ntrks]/I")
outputTree.Branch("tPDG", _tPDG, "tPDG[ntrks]/I")
outputTree.Branch("tEvt", _tEvt, "tEvt[ntrks]/I")
outputTree.Branch("tMot", _tMot, "tMot[ntrks]/I")
outputTree.Branch("prob", _prob, "prob/F")
outputTree.Branch("maxaperture", _maxaperture, "maxaperture/F")
outputTree.Branch("meanaperture", _meanaperture, "meanaperture/F")
outputTree.Branch("maxdphi", _maxdphi, "maxdphi/F")
outputTree.Branch("meanphi", _meanphi, "meanphi/F")
outputTree.Branch("trPDG", _trPDG, "trPDG/F")
outputTree.Branch("signal", _signal, "signal/I")
outputTree.Branch("weight", _weight, "weight/F")
outputTree.Branch("magdphi", _magdphi, "magdphi/F")
outputTree.Branch("dphi", _dphi, "dphi[ntrks]/F")
outputTree.Branch("plate", _plate, "plate/I")
outputTree.Branch("nfirst", _nfirst, "nfirst[ntrks]/I")
outputTree.Branch("nlast", _nlast, "nlast[ntrks]/I")
outputTree.Branch("tLastX", _tLastX, "tLastX[ntrks]/F")
outputTree.Branch("tLastY", _tLastY, "tLastY[ntrks]/F")
outputTree.Branch("tLastZ", _tLastZ, "tLastZ[ntrks]/F")

### VERTICES LOOP ###
vertices = vertextot.eVTX
for ivtx, vtx in enumerate(vertices):
    ntrks = vtx.N()
    vx = vtx.VX()
    vy = vtx.VY()
    vz = vtx.VZ()
    flag = vtx.Flag()
    ntrks = vtx.N()

    apeList = []
    phiList = []
    TXList = []
    TYList = []
    TXNList = []
    TYNList = []
    ffList = []
    ipList = []
    plateList = []

    ### TRACKS LOOP ###
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        nseg = track.N()
        ### SEGMENTS LOOP ###
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
        
        npl = track.Npl()
        tx = track.TX()
        ty = track.TY()
        sx = track.GetSegmentLast().X()
        sy = track.GetSegmentLast().Y()
        sz = track.GetSegmentLast().Z()
        nfirst = track.GetSegmentFirst().Plate()
        nlast = track.GetSegmentLast().Plate()
        plateList.append(nfirst)
        impact_parameter = vtx.GetVTa(itrack).Imp()
        ipList.append(impact_parameter)
        nfirst = track.GetSegmentFirst().Plate()
        nava = from_plate - nfirst + 1
        FF = float(nseg)/float(nava)
        phi = track.Phi()
        ffList.append(FF)
        TNorm = ROOT.TMath.Sqrt(tx*tx + ty*ty)
        TXList.append(tx)
        TYList.append(ty)
        TXNList.append(tx/TNorm)
        TYNList.append(ty/TNorm)
        phiList.append(phi)
        _nseg[itrack] = nseg
        _npl[itrack] = npl
        _fillfact_t[itrack] = FF
        _ip[itrack] = impact_parameter
        _tTX[itrack] = tx
        _tTY[itrack] = ty
        _tPDG[itrack] = trackPDG
        _tEvt[itrack] = trackEvt
        _tID[itrack] = trackID
        _tMot[itrack] = trackMother
        _nfirst[itrack] = nfirst
        _nlast[itrack] = nlast
        _tLastX[itrack] = sx
        _tLastY[itrack] = sy           
        _tLastZ[itrack] = sz   
        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(ROOT.TMath.Sqrt( tx*tx+ty*ty ))
    
    plate = min(plateList)
    _plate[0] = plate

    arrPhi = np.array(phiList)
    arrTX = np.array(TXList)
    arrTY = np.array(TYList)
    arrTXN = np.array(TXNList)
    arrTYN = np.array(TYNList)
    dPhiList = []
    for itrack in range(ntrks):
        average_of_rest = ROOT.TMath.ATan2(np.sum(arrTYN) - arrTYN[itrack], np.sum(arrTXN) - arrTXN[itrack])
        difference = ROOT.TMath.Abs(arrPhi[itrack] - average_of_rest)
        if difference > ROOT.TMath.Pi():
            difference = 2*ROOT.TMath.Pi() - difference
        dPhiList.append(difference)
        _dphi[itrack]=difference
    magdphi = ROOT.TMath.Sqrt(np.sum(arrTX)**2 + np.sum(arrTY)**2)

    _brickID[0] = brickID
    _npart[0] = npart
    _vx[0]=vx
    _vy[0]=vy
    _vz[0]=vz
    _flag[0]=flag
    _vID[0] = vtx.ID()
    _ntrks[0] = ntrks
    _fillfact[0] = np.mean(ffList)
    _meanIP[0] = np.mean(ipList)
    _prob[0] = vtx.V().prob()
    _maxaperture[0] = vtx.MaxAperture()
    _maxdphi[0] = np.max(dPhiList)
    _meanphi[0] = np.mean(phiList)
    _meanaperture[0] = np.mean(apeList)
    _signal[0] = 0
    _weight[0] = w
    outputTree.Fill()
    
del gAli

#write output files
outputFile.cd()
outputTree.Write()
outputFile.Close()

print(f'Reco vtx: {rec_vtx}')
print(f'Nu vtx: {neu_vtx}')
print(f'Nu vtx: {ivtx+1}')

elapsed_time = time() - start_time
print("TOTAL ELAPSED TIME ", elapsed_time)
print("Weight used: ", w)