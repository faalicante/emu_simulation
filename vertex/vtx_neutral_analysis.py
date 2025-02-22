import ROOT
import numpy as np
import fedrarootlogon
from time import time
import os.path
from array import array


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-p", dest="ProcId", required=True, default=None, type=int)
options = parser.parse_args()

nbrick = options.ProcId%5
npart = options.ProcId//5
brickIDs = [11,21,31,41,51]
brickID = brickIDs[nbrick]


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


def prepareList(cbmsim):
    neu_list = []
    for ievt, event in enumerate(cbmsim):
        if ievt%100==0: print(f"Event {ievt}")
        neu_vtx = ROOT.TVector3(event.MCTrack[1].GetStartX(), event.MCTrack[1].GetStartY(), event.MCTrack[1].GetStartZ())
        neu_list.append(neu_vtx)
    return neu_list

def findVertex(vtx, nu_list):
    max_xy=0.0003 #3um
    max_z=0.1350 #1350um
    vtx_conv = convertVertex(vtx, brickID)
    vtx_g = ROOT.TVector3(vtx_conv[0], vtx_conv[1], vtx_conv[2])
    closestEvent = None
    for ievt, nu_vtx in enumerate(nu_list):
        dist_xy = ROOT.TMath.Sqrt((nu_vtx.X()-vtx_g.X())**2 + (nu_vtx.Y()-vtx_g.Y())**2)
        dist_z = ROOT.TMath.Abs(nu_vtx.Z()-vtx_g.Z())
        if dist_z <= max_z:
            if dist_xy < max_xy:
                max_xy = dist_xy
                closestEvent = ievt
    return closestEvent

nbrick = options.ProcId%20
npart = options.ProcId//20
brickIDs = [11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44,51,52,53,54]
brickID = brickIDs[nbrick]
part = ['K_10_50', 'K_50_100', 'K_100_200', 'neu_10_50', 'neu_50_100', 'neu_100_200']

##WEIGHT USED FOR DATA IN RUN1 BRICK21
events_b21 = [181, 206, 203, 195, 189, 214]
expected_b21 = [42.840602, 5.7650435, 2.4583636, 24.623528, 2.1840258, 0.40455220]
lumi = 9.512
fid_area = 0.6917317708
w = expected_b21[npart]*lumi*fid_area/events_b21[npart]

from_plate = 60
zmin = -77585.00

emureader = ROOT.EmulsionDet()


path = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/neutrals/'+part[npart]+'/Ntuples'
out_dir = path
out_name = f'/vertex_neutral_{npart}_{brickID}.root'

prepath = '/eos/experiment/sndlhc/MonteCarlo/FEDRA/neutrals/'
energyBins = ['_10_50', '_50_100', '_100_200']

if npart < 3:
  simFile = path+'/sndLHC.PG_130-TGeant4.root'
  geoFile =  path+'/geofile_full.PG_130-TGeant4.root'
else:
  simFile = path+'/sndLHC.PG_2112-TGeant4.root'
  geoFile =  path+'/geofile_full.PG_2112-TGeant4.root'

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(geoFile)

fsim = ROOT.TFile.Open(simFile)
cbmsim = fsim.cbmsim
neu_list = prepareList(cbmsim)


out_name = f'/vertex_neutral_{npart}_{brickID}.root'
vtx_file = path+'/b{:06d}/b{:06d}.0.0.0.vtx.root'.format(brickID, brickID)


#histo setup
h_n = ROOT.TH1D('n','Multiplicity;multiplicity', 30, 0, 30)
h_flag = ROOT.TH1D('flag','Flag;flag', 6, 0, 6)
h_vz = ROOT.TH1D('vz','Vertex z position;vz[um]', 400, -80000, 5000)
h_vxy = ROOT.TH2D('vxy', 'Vertex xy map;vx[um];vy[um]', 200, 0, 200000, 200, 0, 200000)
h_n0 = ROOT.TH1D('n0', 'Multiplicity;multiplicity', 27, 3, 30)
h_nseg = ROOT.TH1D('nseg', 'Number of segments;nseg', 56, 4, 60)
h_npl = ROOT.TH1D('npl', ' Number of crossing films;npl', 56, 4, 60)
h_ff = ROOT.TH1D('ff', 'Fill Factor;FF', 55, 0, 1.1)
h_ip = ROOT.TH1D('ip', 'Impact parameter;ip[um]', 100, 0, 20)
h_meanff = ROOT.TH1D('meanff', 'Mean Fill Factor;FF', 22, 0, 1.05)
h_meanip = ROOT.TH1D('meanip', 'Mean impact parameter;ip[um]', 100, 0, 20)
h_prob = ROOT.TH1D('prob', 'Probability;prob', 55, 0, 1.1)
h_maxape = ROOT.TH1D('maxape', 'Max aperture;max_ape', 100, 0, 1)
h_meanape = ROOT.TH1D('meanape', 'Mean aperture;mean_ape', 100, 0, 1)
h_meanphi = ROOT.TH1D('meanphi', 'Mean phi;mean_phi', 160, -4, 4)
h_maxdphi = ROOT.TH1D('maxdphi', 'Max phi diff;max_dphi', 80, 0, 4)
h_offset_xy = ROOT.TH2D('offset_xy', 'True neutrino vs vtx;x;y', 200, -0.001, 0.001, 200, -0.001, 0.001)
h_offset_z = ROOT.TH1D('offset_z', 'True neutrino vs vtx;z', 2000, -0.1, 0.1)

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
 
fake_vtx=0

### VERTICES LOOP ###
for ivtx, vtx in enumerate(vertices):
    ntrks = vtx.N()
    neu_vtx=0
    lep_found=0
    fake_tracks=0
    vx = vtx.VX()
    vy = vtx.VY()
    vz = vtx.VZ()
    flag = vtx.Flag()
    ntrks = vtx.N()
    h_vxy.Fill(vx, vy)
    if vx < 0 or vx > 195000: continue
    if vy < 0 or vy > 195000: continue
    h_vz.Fill(vz)
    if vz < zmin or vz > 0: continue
    h_flag.Fill(flag)
    if flag !=0 and flag !=3: continue
    # print(f"Vertex {ivtx}")
    h_n.Fill(ntrks, w)
    if ntrks < 3: continue
    closestEvent = findVertex(vtx, neu_list)
    if closestEvent == None:
        print("no closest event found")
        continue

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
        if trackMother == 0: neu_vtx += 1
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
        DictTrackPdg[track.MCTrack()] = trackPDG
        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(ROOT.TMath.Sqrt( tx*tx+ty*ty ))

    if neu_vtx<0.5*ntrks:
        print("no mother id ",  vtx.ID())
        fake_vtx+=1
        continue

    eventID = max(DictTrackEvt, key=DictTrackEvt.get)

    h_n0.Fill(ntrks, w)
    h_prob.Fill(vtx.V().prob(), w)
    h_maxape.Fill(vtx.MaxAperture(), w)
    h_meanff.Fill(np.mean(ffList), w)
    h_meanip.Fill(np.mean(ipList), w)
    h_meanape.Fill(np.mean(apeList), w)
    h_meanphi.Fill(np.mean(phiList), w)

    nx = cbmsim.MCTrack[1].GetStartX()
    ny = cbmsim.MCTrack[1].GetStartY()
    nz = cbmsim.MCTrack[1].GetStartZ()

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
    _npart[0] = npart
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
    _signal[0] = 0
    _weight[0] = w
    outputTree.Fill()
    
del gAli

#write output files
outputFile.cd()
outputTree.Write()
outputFile.Close()

histoFile = ROOT.TFile(out_dir+f"/hist_out_{npart}_{brickID}.root", "RECREATE")
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
elapsed_time = time() - start_time
print("TOTAL ELAPSED TIME ", elapsed_time)
print("Weight used: ", w)