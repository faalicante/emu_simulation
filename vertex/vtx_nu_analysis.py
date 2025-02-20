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
parser.add_argument("--muon", dest="muon", required=False, default=None)
options = parser.parse_args()

numu = options.numu
nue = options.nue
muon = options.muon
brickID = options.brickID
from_plate = 60

def AddToDict(dictionary, var):
    if var in dictionary:
        dictionary[var] += 1
    else:
        dictionary[var] = 1

def evalDiffPhi(vtx, phiList, TXList, TYList):
    arrPhi = np.array(phiList)
    arrTX = np.array(TXList)
    arrTY = np.array(TYList)
    dPhiList = {}
    for itrack in range(vtx.N()):
        track = vtx.GetTrack(itrack)
        trackID = track.MCTrack()
        aor = ROOT.TMath.ATan2(np.sum(arrTY) - arrTY[itrack], np.sum(arrTX) - arrTX[itrack])
        difference = ROOT.TMath.Abs(arrPhi[itrack] - aor)
        if difference > ROOT.TMath.Pi():
            difference = 2*ROOT.TMath.Pi() - difference
        dPhiList[trackID] = difference
    return dPhiList

def convertVertex(vtx,brickID,refplate = 60):
    # convert position and angles from FEDRA into the global system
    detID = int(brickID * 1e3 + refplate)
    localarr = np.array([vtx.VX(), vtx.VY(), vtx.VZ()])
    globalarr = np.array([0., 0., 0.])
    emureader.GetPosition(detID,localarr,globalarr)
    return globalarr

def checkDict(dictionary):
    miss = 0
    for key, value in dictionary.items():
        if key != value: miss+=1
    return miss

def prepareList(cbmsim):
    nu_list = []
    for ievt, event in enumerate(cbmsim):
        if ievt%100==0: print(f"Event {ievt}")
        nu_vtx = ROOT.TVector3(event.MCTrack[0].GetStartX(), event.MCTrack[0].GetStartY(), event.MCTrack[0].GetStartZ())
        nu_list.append(nu_vtx)
    return nu_list

def findVertex(vtx, nu_list):
    max_xy=0.0003 #3um
    max_z=0.1350 #1350um
    vtx_conv = convertVertex(vtx, brickID)
    vtx_g = ROOT.TVector3(vtx_conv[0], vtx_conv[1], vtx_conv[2])
    closestEvent = None
    for ievt, nu_vtx in enumerate(nu_list):
        dist_xy = ROOT.TMath.Sqrt((nu_vtx.X()-vtx_g.X())**2 + (nu_vtx.Y()-vtx_g.Y())**2)
        dist_z = ROOT.TMath.Abs(nu_vtx.Z()-vtx_g.Z())
        if dist_z < max_z:
            # max_z = dist_z
            if dist_xy < max_xy:
                # max_xy = dist_xy
                closestEvent = ievt
    return closestEvent


zmin = -77585.00
# pathSim = '/eos/experiment/sndlhc/MonteCarlo/Neutrinos/Genie/nu_sim_activeemu_withcrisfiles_25_July_2022/'

if numu:
    # path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc_eff9_smear1'
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc_eff10_smear0'
    sim_file = path+'/inECC_sndLHC.Genie-TGeant4.root'
    out_name = f'/vertex_sigmu_dan_{brickID}.root'
elif nue:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc_eff9_smear1'
    # path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/nuecc_eff10_smear0'
    sim_file = path+'/inECC_sndLHC.Genie-TGeant4.root'
    out_name = f'/vertex_sige_{brickID}.root'
elif muon:
    path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1E5_eff9_smear1'
    sim_file = path+'/sndLHC.Ntuple-TGeant4-1E5cm2.root'
    out_name = f'/vertex_muon_{brickID}.root'

geoFile =  path + '/geofile_full.Genie-TGeant4.root'
import SndlhcGeo
geo = SndlhcGeo.GeoInterface(geoFile)

out_dir = path
vtx_file = path+'/b{:06d}/b{:06d}.0.0.0.vtx.root'.format(brickID, brickID)

#histo setup
h_n = ROOT.TH1D('n','Multiplicity;multiplicity', 50, 0, 50)
h_flag = ROOT.TH1D('flag','Flag;flag', 6, 0, 6)
h_vz = ROOT.TH1D('vz','Vertex z position;vz[um]', 400, -80000, 5000)
h_vxy = ROOT.TH2D('vxy', 'Vertex xy map;vx[um];vy[um]', 200, 0, 200000, 200, 0, 200000)
h_n0 = ROOT.TH1D('n0', 'Multiplicity;multiplicity', 50, 0, 50)
h_nseg = ROOT.TH1D('nseg', 'Number of segments;nseg', 56, 4, 60)
h_npl = ROOT.TH1D('npl', ' Number of crossing films;npl', 56, 4, 60)
h_ff = ROOT.TH1D('ff', 'Fill Factor;FF', 22, 0, 1.05)
h_ip = ROOT.TH1D('ip', 'Impact parameter;ip[um]', 500, 0, 50)
h_meanff = ROOT.TH1D('meanff', 'Mean Fill Factor;FF', 22, 0, 1.05)
h_meanip = ROOT.TH1D('meanip', 'Mean impact parameter;ip[um]', 500, 0, 50)
h_prob = ROOT.TH1D('prob', 'Probability;prob', 30, 0, 1.02)
h_maxape = ROOT.TH1D('maxape', 'Max aperture;max_ape', 50, 0, 2.5)
h_meanape = ROOT.TH1D('meanape', 'Mean aperture;mean_ape', 50, 0, 2.5)
h_meanphi = ROOT.TH1D('meanphi', 'Mean phi;mean_phi', 80, -4, 4)
h_maxdphi = ROOT.TH1D('maxdphi', 'Max phi diff;max_dphi', 40, 0, 4)
h_offset_xy = ROOT.TH2D('offset_xy', 'True neutrino vs vtx;x;y', 200, -0.001, 0.001, 200, -0.001, 0.001)
h_offset_z = ROOT.TH1D('offset_z', 'True neutrino vs vtx;z', 2000, -0.1, 0.1)

emureader = ROOT.EmulsionDet()

N=40
start_time = time()
#reading vertices
if not os.path.isfile(vtx_file):
    print(f"{vtx_file} not found, interrupting")
    exit()
#save vertices in root file
outputFile = ROOT.TFile(out_dir+out_name,"RECREATE")	
outputTree = ROOT.TTree("vtx","Tree of vertices")
outputTree.SetDirectory(outputFile)

_brickID = array('i', [0])
_vID = array('i', [0])
_flag = array('i', [0])
_vx = array('f', [0])
_vy = array('f', [0])
_vz = array('f', [0])
_vx_g = array('f', [0])
_vy_g = array('f', [0])
_vz_g = array('f', [0])
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
_maxdphi = array('f', [0])
_meanphi = array('f', [0])
_tTX = array('f', N*[0])
_tTY = array('f', N*[0])
_tID = array('i', N*[0])
_tPDG = array('i', N*[0])
_tEvt = array('i', N*[0])
_MCEvt = array('i', [0])
_trPDG = array('f', [0])
_motherdPhi = array('f', [0])
_signal = array('i', [0])
_weight = array('f', [0])
_f_trk = array('i', [0])
_clEvt = array('i', [0])

outputTree.Branch("brickID", _brickID, "brickID/I")
outputTree.Branch("vID", _vID, "vID/I")
outputTree.Branch("flag", _flag, "flag/I")
outputTree.Branch("vx_g", _vx_g, "vx_g/F")
outputTree.Branch("vy_g", _vy_g, "vy_g/F")
outputTree.Branch("vz_g", _vz_g, "vz_g/F")
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
outputTree.Branch("prob", _prob, "prob/F")
outputTree.Branch("maxaperture", _maxaperture, "maxaperture/F")
outputTree.Branch("meanaperture", _meanaperture, "meanaperture/F")
outputTree.Branch("maxdphi", _maxdphi, "maxdphi/F")
outputTree.Branch("meanphi", _meanphi, "meanphi/F")
outputTree.Branch("MCEvt", _MCEvt, "MCEvt/I")
outputTree.Branch("trPDG", _trPDG, "trPDG/F")
outputTree.Branch("motherdPhi", _motherdPhi, "motherdPhi/F")
outputTree.Branch("signal", _signal, "signal/I")
outputTree.Branch("weight", _weight, "weight/F")
outputTree.Branch("f_trk", _f_trk, "f_trk/I")
outputTree.Branch("clEvt", _clEvt, "clEvt/I")

print("opening file: ",vtx_file)
dproc = ROOT.EdbDataProc()
gAli = dproc.PVR()
scancond = ROOT.EdbScanCond()
scancond.SetSigma0(1.5,1.5,0.0015,0.0015)
scancond.SetDegrad(3)
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
dproc.ReadVertexTree(vertexrec, vtx_file, "1")
vertices = gAli.eVTX

fsim = ROOT.TFile.Open(sim_file)
cbmsim = fsim.cbmsim
nu_list = prepareList(cbmsim)

fake_vtx=0
miss_lep=0
miss_id=0
miss_evt=0
miss_evt2=0
miss_evt_close=0
miss_mother=0
miss_pdg=0

### VERTICES LOOP ###
for ivtx, vtx in enumerate(vertices):
    ntrks = vtx.N()
    nu_vtx=0
    lep_found=0
    fake_tracks=0
    vx = vtx.VX()
    vy = vtx.VY()
    vz = vtx.VZ()
    flag = vtx.Flag()
    ntrks = vtx.N()
    h_vz.Fill(vz)
    h_vxy.Fill(vx, vy)
    if vz < zmin or vz > 0: continue
    if vx < 0 or vx > 200000: continue
    if vy < 0 or vy > 200000: continue
    h_flag.Fill(flag)
    if flag !=0 and flag !=3: continue
    # print(f"Vertex {ivtx}")
    closestEvent = findVertex(vtx, nu_list)
    if closestEvent == None:
        # fake_vtx+=1
        continue
    h_n.Fill(ntrks)
    # if ntrks < 3: continue

    apeList = []
    phiList = []
    TXList = []
    TYList = []
    ffList = []
    ipList = []

    ### TRACKS LOOP ###
    DictTrackPdg = {}
    DictTrackEvt = {}
    DictCheckEvt = {}
    DictCheckEvtClose = {}
    DictCheckID = {}
    DictCheckMother = {}
    DictCheckPdg = {}
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
            AddToDict(DictSegEvt, sEvt)
            AddToDict(DictSegID, sID)
            AddToDict(DictSegMother, sMother)
        trackPDG = max(DictSegPDG, key=DictSegPDG.get)
        trackEvt = max(DictSegEvt, key=DictSegEvt.get)
        trackID = max(DictSegID, key=DictSegID.get)
        trackMother = max(DictSegMother, key=DictSegMother.get)

        # DictCheckID[trackID] = track.MCTrack()
        # DictCheckID2[trackID] = track.GetSegmentsMCTrack(nseg_int)
        # DictCheckEvt[trackEvt] = track.MCEvt()
        # DictCheckMother[trackMother] = track.Aid(0)
        # DictCheckPdg[trackPDG] = track.GetSegment(0).Vid(0)
        # track_out = vtx.GetVTa(itrack).Zpos()
        if trackMother == 0:
            if trackID != track.MCTrack():
                print("id", vtx.ID(), trackID, track.MCTrack())
                miss_id+=1
            if trackEvt != track.MCEvt():
                print("evt", vtx.ID(), trackID, track.MCTrack())
                miss_evt+=1
            if trackMother != track.Aid(0):
                print("mot", vtx.ID(), trackID, track.MCTrack())
                miss_mother+=1
            if trackPDG != track.GetSegment(0).Vid(0):
                print("pdg", vtx.ID(), trackID, track.MCTrack())
                miss_pdg+=1
            nu_vtx += 1
            if numu and abs(trackPDG) == 13:
                lep_found = True
            if nue and abs(trackPDG) == 11:
                lep_found = True
        #     # nu_vtx += 1
        #     elif nue and abs(track.GetSegment(0).Vid(0)) == 11:
        #         nu_vtx = True
        npl = track.Npl()
        impact_parameter = vtx.GetVTa(itrack).Imp()
        h_ip.Fill(impact_parameter)         
        ipList.append(impact_parameter)
        nfirst = track.GetSegmentFirst().Plate()
        nava = from_plate - nfirst + 1
        FF = float(nseg)/float(nava)
        ffList.append(FF)
        TXList.append(track.TX())
        TYList.append(track.TY())
        phiList.append(track.Phi())
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
        DictTrackPdg[track.MCTrack()] = trackPDG
        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(ROOT.TMath.Sqrt( tx*tx+ty*ty ))

    # miss_id += checkDict(DictCheckID)
    # miss_id2 += checkDict(DictCheckID2)
    # miss_evt += checkDict(DictCheckEvt)
    # miss_mother += checkDict(DictCheckMother)
    # miss_pdg += checkDict(DictCheckPdg)

    # if nu_vtx<0.5*ntrks: 
    if (numu or nue) and nu_vtx<0.5*ntrks:
        print("no mother id ",  vtx.ID())
        fake_vtx+=1
        continue
    if not lep_found:
        print("no lep" , vtx.ID())
        miss_lep+=1
    fake_tracks = ntrks - nu_vtx
    _f_trk[0] = fake_tracks
    eventID = max(DictTrackEvt, key=DictTrackEvt.get)
    if eventID != closestEvent:
        print("cl ev", vtx.ID(), eventID, closestEvent)
        miss_evt_close +=1
    # DictCheckEvtClose[eventID] = closestEvent
    # miss_evt_close += checkDict(DictCheckEvtClose)
    cbmsim.GetEntry(eventID)
    # eventID = cbmsim.MCEventHeader.GetEventID()
    w = cbmsim.MCTrack[0].GetWeight()
    h_n0.Fill(ntrks)
    h_prob.Fill(vtx.V().prob())
    h_maxape.Fill(vtx.MaxAperture())
    h_meanff.Fill(np.mean(ffList))
    h_meanip.Fill(np.mean(ipList))
    h_meanape.Fill(np.mean(apeList))
    h_meanphi.Fill(np.mean(phiList))

    nx = cbmsim.MCTrack[0].GetStartX()
    ny = cbmsim.MCTrack[0].GetStartY()
    nz = cbmsim.MCTrack[0].GetStartZ()

    vtx_g = convertVertex(vtx, brickID)
    vx_g = vtx_g[0]
    vy_g = vtx_g[1]
    vz_g = vtx_g[2]
    h_offset_xy.Fill(nx-vx_g, ny-vy_g)
    h_offset_z.Fill(nz-vz_g)

    dPhiList = evalDiffPhi(vtx, phiList, TXList, TYList)
    track_maxdphi = max(dPhiList, key=dPhiList.get)
    maxdphi = dPhiList[track_maxdphi]
    h_maxdphi.Fill(maxdphi)

    _brickID[0] = brickID
    _trPDG[0] = DictTrackPdg[track_maxdphi]
    _MCEvt[0] = eventID
    _clEvt[0] = closestEvent
    if track_maxdphi < len(cbmsim.MCTrack):
        motherID = cbmsim.MCTrack[track_maxdphi].GetMotherId()
        _motherdPhi[0] = motherID
    else:
        _motherdPhi[0] = 10
    
    _vx[0]=vx
    _vy[0]=vy
    _vz[0]=vz
    _vx_g[0]=vx_g
    _vy_g[0]=vy_g
    _vz_g[0]=vz_g
    _flag[0]=flag
    _vID[0] = vtx.ID()
    _ntrks[0] = ntrks
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
outputFile.cd()
outputTree.Write()
outputFile.Close()

histoFile = ROOT.TFile(out_dir+"/hist_out_dan_{}.root".format(brickID), "RECREATE")
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
h_offset_xy.Write()
h_offset_z.Write()
histoFile.Write()
histoFile.Close()

print(f"fake vertices: {fake_vtx}")
print(f"missing lepton: {miss_lep}")
print("Missing track ID: ", miss_id)
print("Missing track Mother: ", miss_mother)
print("Missing track PDG: ", miss_pdg)
print("Missing track Evt: ", miss_evt)
print("Missing track Evt Close: ", miss_evt_close)
elapsed_time = time() - start_time
print("TOTAL ELAPSED TIME ", elapsed_time)