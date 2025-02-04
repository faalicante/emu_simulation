import os.path
from time import time
from array import array
import ROOT as r
import numpy as np
import fedrarootlogon

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-b", dest="brickID", required=True, default=None, type=int)
parser.add_argument("--numu", dest="numu", required=False, default=None)
parser.add_argument("--nue", dest="nue", required=False, default=None)
parser.add_argument("--muon", dest="muon", required=False, default=None)
options = parser.parse_args()

numu = options.numu
nue = options.nue
muon = options.muon
brickID = options.brickID
from_plate = 60

zmin = -77585.00

if numu:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc_eff9_smear1'
    sim_file = path+'/inECC_sndLHC.Genie-TGeant4.root'
    out_name = f'/vertex_sigmu_{brickID}.root'
elif nue:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc_eff9_smear1'
    sim_file = path+'/inECC_sndLHC.Genie-TGeant4.root'
    out_name = f'/vertex_sige_{brickID}.root'
elif muon:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1E5_eff9_smear1'
    sim_file = path+'/sndLHC.Ntuple-TGeant4-1E5cm2.root'
    out_name = f'/vertex_muon_{brickID}.root'
out_dir = path
vtx_file = path+'/b{:06d}/b{:06d}.0.0.0.vtx.root'.format(brickID, brickID)

#histo setup
h_n = r.TH1D('n','Multiplicity;multiplicity', 50, 0, 50)
h_flag = r.TH1D('flag','Flag;flag', 6, 0, 6)
h_vz = r.TH1D('vz','Vertex z position;vz[um]', 400, -80000, 5000)
h_vxy = r.TH2D('vxy', 'Vertex xy map;vx[um];vy[um]', 200, 0, 200000, 200, 0, 200000)
h_n0 = r.TH1D('n0', 'Multiplicity;multiplicity', 47, 3, 50)
h_nseg = r.TH1D('nseg', 'Number of segments;nseg', 56, 4, 60)
h_npl = r.TH1D('npl', ' Number of crossing films;npl', 56, 4, 60)
h_ff = r.TH1D('ff', 'Fill Factor;FF', 22, 0, 1.05)
h_ip = r.TH1D('ip', 'Impact parameter;ip[um]', 100, 0, 300)
h_meanff = r.TH1D('meanff', 'Mean Fill Factor;FF', 22, 0, 1.05)
h_meanip = r.TH1D('meanip', 'Mean impact parameter;ip[um]', 100, 0, 300)
h_prob = r.TH1D('prob', 'Probability;prob', 30, 0, 1.02)
h_maxape = r.TH1D('maxape', 'Max aperture;max_ape', 50, 0, 2.5)
h_meanape = r.TH1D('meanape', 'Mean aperture;mean_ape', 50, 0, 2.5)
h_meanphi = r.TH1D('meanphi', 'Mean phi;mean_phi', 80, -4, 4)
h_maxdphi = r.TH1D('maxdphi', 'Max phi diff;max_dphi', 40, 0, 4)
h_offset = r.TH2D('offset', 'True neutrino vs vtx;x;y', 200, -20000, 20000, 200, -20000, 20000)

N=40
start_time = time()
#reading vertices
if not os.path.isfile(vtx_file):
    print(f"{vtx_file} not found, interrupting")
    exit()
#save vertices in root file
outputTree = r.TTree("vtx","Tree of vertices")

_brickID = array('i', [0])
_vID = array('i', [0])
_flag = array('i', [0])
_vx = array('f', [0])
_vy = array('f', [0])
_vz = array('f', [0])
_ntrks = array('i', [0])
_nsegtot = array('i', [0])
_nseg = array('i', N*[0])
_npl = array('i', N*[0])
_fillfact_t = array('f', N*[0])
_fillfact = array('f', [0])
_meanIP = array('f', [0])
_ip = array('f', N*[0])
_prob = array('f', [0])
_maxaperture = array('f', [0])
_meanaperture = array('f', [0])
_maxdphi = array('f', [0])
_meanphi = array('f', [0])
_tTX = array('f', N*[0])
_tTY = array('f', N*[0])
_tID = array('i', N*[0])
_tPDG = array('i', N*[0])
_tEvt = array('i', N*[0])
_MCEvt = array('i', [0])
_evtPDG = array('f', [0])
_trPDG = array('f', [0])
_motherdPhi = array('f', [0])
_signal = array('i', [0])
_weight = array('f', [0])

outputTree.Branch("brickID", _brickID, "brickID/I")
outputTree.Branch("vID", _vID, "vID/I")
outputTree.Branch("flag", _flag, "flag/I")
outputTree.Branch("vx", _vx, "vx/F")
outputTree.Branch("vy", _vy, "vy/F")
outputTree.Branch("vz", _vz, "vz/F")
outputTree.Branch("ntrks", _ntrks, "ntrks/I")
outputTree.Branch("nsegtot", _nsegtot, "nsegtot/I")
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
outputTree.Branch("prob", _prob, "prob/F")
outputTree.Branch("maxaperture", _maxaperture, "maxaperture/F")
outputTree.Branch("meanaperture", _meanaperture, "meanaperture/F")
outputTree.Branch("maxdphi", _maxdphi, "maxdphi/F")
outputTree.Branch("meanphi", _meanphi, "meanphi/F")
outputTree.Branch("MCEvt", _MCEvt, "MCEvt/I")
outputTree.Branch("evtPDG", _evtPDG, "evtPDG/F")
outputTree.Branch("trPDG", _trPDG, "trPDG/F")
outputTree.Branch("motherdPhi", _motherdPhi, "motherdPhi/F")
outputTree.Branch("signal", _signal, "signal/I")
outputTree.Branch("weight", _weight, "weight/F")

print("opening file: ",vtx_file)
dproc = r.EdbDataProc()
gAli = dproc.PVR()
scancond = r.EdbScanCond()
scancond.SetSigma0(1.5,1.5,0.0015,0.0015)
scancond.SetDegrad(3)
gAli.SetScanCond(scancond)
vertexrec = r.EdbVertexRec()
vertexrec.SetPVRec(gAli)
vertexrec.eDZmax=3000.
vertexrec.eProbMin=0.0001
vertexrec.eImpMax=15.
vertexrec.eUseMom=False
vertexrec.eUseSegPar=True
vertexrec.eQualityMode=0
proc = r.EdbDataProc()
dproc.ReadVertexTree(vertexrec, vtx_file, "1")
vertices = gAli.eVTX

fsim = r.TFile.Open(sim_file)
cbmsim = fsim.cbmsim

fake_vtx=0
### VERTICES LOOP ###
for vtx in vertices:
    ntrks = vtx.N()
    nu_vtx=0
    vx = vtx.VX()
    vy = vtx.VY()
    vz = vtx.VZ()
    flag = vtx.Flag()
    ntrks = vtx.N()
    h_vz.Fill(vz)
    h_vxy.Fill(vx, vy)
    h_flag.Fill(flag)
    h_n.Fill(ntrks)
    if vz < zmin: continue
    if vz > 0: continue
    if flag !=0 and flag !=3: continue
    if ntrks < 3: continue

    apeList = []
    phiList = []
    TXList = []
    TYList = []
    ffList = []
    ipList = []
    segidx = 0

    ### TRACKS LOOP ###
    DictTrackPdg = {}
    DictTrackEvt = {}
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        tEvt = track.MCEvt()
        tID = track.MCTrack()
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
        track_out = vtx.GetVTa(itrack).Zpos()
        if trackMother == 0:
            if numu and abs(trackPDG) == 13:
                nu_vtx = True
            elif nue and abs(trackPDG) == 11:
                nu_vtx = True
        npl = track.Npl()
        impact_parameter = vtx.GetVTa(itrack).Imp()
        h_ip.Fill(impact_parameter)         
        ipList.append(impact_parameter)
        nfirst = track.GetSegmentFirst().Plate()
        nava = from_plate - nfirst + 1
        FF = float(nseg)/float(nava)
        phi = track.Phi()
        ffList.append(FF)
        TXList.append(track.TX())
        TYList.append(track.TY())
        phiList.append(phi)
        h_nseg.Fill(nseg)
        h_npl.Fill(npl)
        h_ff.Fill(FF)
        _nseg[itrack] = nseg
        _npl[itrack] = npl
        _fillfact_t[itrack] = FF
        _ip[itrack] = impact_parameter
        _tTX[itrack] = track.TX()
        _tTY[itrack] = track.TY()
        _tPDG[itrack] = trackPDG
        _tEvt[itrack] = trackEvt
        _tID[itrack] = trackID
        AddToDict(DictTrackEvt, trackEvt)

        tID = track.MCTrack()
        DictTrackPdg[tID] = trackPDG
        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(r.TMath.Sqrt( tx*tx+ty*ty ))

    if (numu or nue) and not nu_vtx:
        fake_vtx+=1
        continue
    eventID = max(DictTrackEvt, key=DictTrackEvt.get)
    cbmsim.GetEntry(eventID)
    w = cbmsim.MCTrack[0].GetWeight()
    h_n0.Fill(ntrks)
    h_prob.Fill(vtx.V().prob())
    h_maxape.Fill(vtx.MaxAperture())
    h_meanff.Fill(np.mean(ffList))
    h_meanip.Fill(np.mean(ipList))
    h_meanape.Fill(np.mean(apeList))
    h_meanphi.Fill(np.mean(phiList))

    nx = cbmsim.MCTrack[0].GetStartX() * 1E+4 + 273000
    ny = cbmsim.MCTrack[0].GetStartY() * 1E+4 - 158000
    h_offset.Fill(nx-vx, ny-vy)

    arrPhi = np.array(phiList)
    arrTX = np.array(TXList)
    arrTY = np.array(TYList)
    dPhiList = {}
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        trackID = track.MCTrack()
        aor = r.TMath.ATan2(np.sum(arrTY) - arrTY[itrack], np.sum(arrTX) - arrTX[itrack])
        difference = r.TMath.Abs(arrPhi[itrack] - aor)
        if difference > r.TMath.Pi():
            difference = 2*r.TMath.Pi() - difference
        dPhiList[trackID] = difference
    track_maxdphi = max(dPhiList, key=dPhiList.get)
    maxdphi = dPhiList[track_maxdphi]
    h_maxdphi.Fill(maxdphi)

    _brickID[0] = brickID
    _trPDG[0] = DictTrackPdg[track_maxdphi]
    _MCEvt[0] = eventID
    if track_maxdphi < len(cbmsim.MCTrack):
        motherID = cbmsim.MCTrack[track_maxdphi].GetMotherId()
        _motherdPhi[0] = motherID
    else:
        _motherdPhi[0] = 10
    
    _vx[0]=vx
    _vy[0]=vy
    _vz[0]=vz
    _flag[0]=flag
    _vID[0] = vtx.ID()
    _ntrks[0] = ntrks
    _nsegtot[0] = segidx
    _fillfact[0] = np.mean(ffList)
    _meanIP[0] = np.mean(ipList)
    _prob[0] = vtx.V().prob()
    _maxaperture[0] = vtx.MaxAperture()
    _maxdphi[0] = maxdphi
    _meanphi[0] = np.mean(phiList)
    _meanaperture[0] = np.mean(apeList)
    _signal[0] = 1
    _weight[0] = w
    outputTree.Fill()
    
fsim.Close()
del gAli

#write output files
outputFile = r.TFile(out_dir+out_name.format(brickID),"RECREATE")	
outputTree.Write()

histoFile = r.TFile(out_dir+"/hist_out_{}.root".format(brickID), "RECREATE")
h_n.Write()
h_flag.Write()
h_vxy.Write()
h_vz.Write()
h_n0.Write()
h_nseg.Write()
h_npl.Write()
h_ff.Write()
h_ip.Write()
h_meanff.Write()
h_meanip.Write()
h_prob.Write()
h_maxape.Write()
h_meanape.Write()
h_maxdphi.Write()
h_meanphi.Write()
h_offset.Write()
histoFile.Write()
histoFile.Close()
outputFile.Close()

print(f"fake vertices: {fake_vtx}")
elapsed_time = time() - start_time
print("TOTAL ELAPSED TIME ", elapsed_time)