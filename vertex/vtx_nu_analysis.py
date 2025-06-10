import os.path
from time import time
from array import array
import ROOT
import numpy as np
import fedrarootlogon
from ctypes import c_int

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b", dest="brickID", required=True, default=None, type=int)
parser.add_argument("--numu", dest="numu", required=False, default=None)
parser.add_argument("--nue", dest="nue", required=False, default=None)
parser.add_argument("-p", dest="partition", required=False, default=None, type=int)
options = parser.parse_args()

mom_est = ROOT.EdbMomentumEstimator()

numu = options.numu
nue = options.nue
brickID = options.brickID
partition = options.partition
from_plate = 60
zmin = -77585.00
zmax = 0

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

if numu:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_muon1.3E5'
    sim_path = '/eos/experiment/sndlhc/users/dancc/NUSIM/numu_inBrick21'
    out_name = f'/vertex_numu_{brickID}_{partition}.root'
elif nue:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc/nuecc_muon1.3E5'
    sim_path = '/eos/experiment/sndlhc/users/dancc/NUSIM/nue_inBrick21'
    out_name = f'/vertex_nue_{brickID}_{partition}.root'

# geoFile =  path + '/geofile_full.Genie-TGeant4.root'
# import SndlhcGeo
# geo = SndlhcGeo.GeoInterface(geoFile)

sTree = ROOT.TChain("cbmsim")
for i in range(partition*10, (partition+1)*10):
    sTree.Add(sim_path+f'/{i+1}/sndLHC.Genie-TGeant4.root')

h_en = ROOT.TH1F("h_en", "Energy distribution of neutrinos", 100, 0, 5000)
for entry in sTree:
    h_en.Fill(entry.MCTrack[0].GetEnergy())

out_dir = path+'/vtx_analysis'
outputFile = ROOT.TFile(out_dir+out_name,"RECREATE")	
outputTree = ROOT.TTree("vtx","Tree of vertices")
outputTree.SetDirectory(outputFile)
# vtx_file = path+'/b{:06d}/b{:06d}.0.0.0.vtx.root'.format(brickID, brickID)

# emureader = ROOT.EmulsionDet()
mom_est = ROOT.EdbMomentumEstimator()

N=40
start_time = time()
#reading vertices
# if not os.path.isfile(vtx_file):
#     print(f"{vtx_file} not found, interrupting")
#     exit()
#save vertices in root file

# print("opening file: ",vtx_file)
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
vtx_cut = f'n>2&&flag==0&&vz>{zmin}&&vz<{zmax}'

rec_vtx = 0
nu_vtx = 0
for i in range(partition*1000, (partition+1)*1000):
    vtx_file = path+f"/b0000{brickID}/b0000{brickID}.0.0.{i+1}.vtx.root"
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
        nu_vtx+=1

_brickID = array('i', [0])
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
_signal = array('i', [0])
_mom = array('f', [0])
_mom_long = array('f', [0])
_mom_max = array('f', [0])
_mom_cat = array('f', [0])
_mom_t = array('f', N*[0])
_energy = array('f', [0])

outputTree.Branch("brickID", _brickID, "brickID/I")
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
outputTree.Branch("signal", _signal, "signal/I")
outputTree.Branch("magdphi", _magdphi, "magdphi/F")
outputTree.Branch("dphi", _dphi, "dphi[ntrks]/F")
outputTree.Branch("plate", _plate, "plate/I")
outputTree.Branch("nfirst", _nfirst, "nfirst[ntrks]/I")
outputTree.Branch("nlast", _nlast, "nlast[ntrks]/I")
outputTree.Branch("tLastX", _tLastX, "tLastX[ntrks]/F")
outputTree.Branch("tLastY", _tLastY, "tLastY[ntrks]/F")
outputTree.Branch("tLastZ", _tLastZ, "tLastZ[ntrks]/F")
outputTree.Branch("mom", _mom, "mom/F")
outputTree.Branch("mom_long", _mom_long, "mom_long/F")
outputTree.Branch("mom_max", _mom_max, "mom_max/F")
outputTree.Branch("mom_cat", _mom_cat, "mom_cat/F")
outputTree.Branch("mom_t", _mom_t, "mom_t[ntrks]/F")
outputTree.Branch("energy", _energy, "energy/F")


# fsim = ROOT.TFile.Open(sim_file)
# cbmsim = fsim.cbmsim
# nu_list = prepareList(cbmsim)
vertices = vertextot.eVTX

##START ANALYSIS
### VERTICES LOOP ###
for ivtx, vtx in enumerate(vertices):
    DictVtxEvt = {}
    vtx_mom = 0
    mom_max = 0
    mom_long = 0
    seg_max = 0
    valid_rec = 0
    cat = 0

    ntrks = vtx.N()
    vx = vtx.VX()
    vy = vtx.VY()
    vz = vtx.VZ()
    flag = vtx.Flag()
    ntrks = vtx.N()
    # closestEvent = findVertex(vtx, nu_list)
    # if closestEvent == None:
    #     # fake_vtx+=1
    #     continue

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
        # tEvt = track.MCEvt()
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
            AddToDict(DictVtxEvt, sEvt)
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
        ffList.append(FF)
        TNorm = ROOT.TMath.Sqrt(tx*tx + ty*ty)
        TXList.append(tx)
        TYList.append(ty)
        TXNList.append(tx/TNorm)
        TYNList.append(ty/TNorm)
        phiList.append(track.Phi())
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

        if nseg < 6: track_mom = -6
        else: track_mom = mom_est.PMScoordinate(track)
        _mom_t[itrack] = track_mom
        if track_mom > 0 and track_mom < 1e9:
            valid_rec += 1
            vtx_mom += track_mom 
            if track_mom > mom_max:
                mom_max = track_mom
            if nseg > seg_max:
                mom_long = track_mom
                seg_max = nseg

        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(ROOT.TMath.Sqrt( tx*tx+ty*ty ))

    fract_valid = float(valid_rec)/float(ntrks)
    if fract_valid <= 0.2:
        cat = 1
    elif fract_valid <= 0.5:
        cat = 2  
    elif fract_valid <= 0.8:
        cat = 3
    else:
        cat = 4

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

    vtxEvt = max(DictVtxEvt, key=DictVtxEvt.get)
    sTree.GetEntry(vtxEvt)
    _energy[0] = sTree.MCTrack[0].GetEnergy()

    _brickID[0] = brickID
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
    _magdphi[0] = magdphi
    _meanphi[0] = np.mean(phiList)
    _meanaperture[0] = np.mean(apeList)
    _signal[0] = 1
    _mom[0] = vtx_mom
    _mom_long[0] = mom_long
    _mom_max[0] = mom_max
    _mom_cat[0] = cat
    outputTree.Fill()
    
del gAli

#write output files
outputFile.cd()
outputTree.Write()
h_en.Write()
outputFile.Close()

print(f'Reco vtx: {rec_vtx}')
print(f'Nu vtx: {nu_vtx}')
print(f'Nu vtx: {ivtx+1}')

elapsed_time = time() - start_time
print("TOTAL ELAPSED TIME ", elapsed_time)