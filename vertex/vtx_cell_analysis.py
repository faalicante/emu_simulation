import ROOT
import numpy as np
import fedrarootlogon
from time import time
import os.path
from array import array
import bisect


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cell", dest="cellID", required=True, default=None, type=int)
parser.add_argument("-b", dest="brick", required=True, default=None, type=int)
options = parser.parse_args()

def right_index(sorted_list, target):
    pos = bisect.bisect_left(sorted_list, target)
    if pos == len(sorted_list):
        return len(sorted_list) - 1
    return pos+1

mom_est = ROOT.EdbMomentumEstimator()

brick = options.brick
ncells = 18
mincell = 5000
maxcell = 185000
from_plate = 57
cellsize= (maxcell-mincell)/ncells

icell = options.cellID
xbin=((icell % ncells)) + 1
ybin=((icell // ncells)) + 1
xcenter = xbin * cellsize
ycenter = ybin * cellsize

# path = '/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN2/b000331/cells'
path = '/eos/experiment/sndlhc/users/dancc/MuonDIS/2024/EmuActive_z_2.88_3.55m_Ioni_latelateFLUKA_smeared/RUN1B21_mirror/b000121/cells'
out_dir = '/afs/cern.ch/work/f/falicant/public/RUN1/b121_edipoda'


# with open(out_dir+"/vertices_of_interest.txt", "r") as file:
#     lines = file.readlines()[1:]
# voiid = []
# for line in lines:
#     values = line.split()
#     if int(values[0])!=xbin*10 or int(values[1])!=ybin*10: continue
#     voiid.append(int(values[2]))

N=50
start_time = time()
cell_path = path+f'/cell_{xbin*10}_{ybin*10}/b{brick:06d}'
vtx_file = cell_path+f'/b{brick:06d}.0.0.0.vtx.refit.root'
#reading vertices
if not os.path.isfile(vtx_file):
    print("Vertex file not found, ", vtx_file)
    exit(1)

print("opening file: ",vtx_file)
#scanset
sspath = cell_path+"/.."
sproc = ROOT.EdbScanProc()
sproc.eProcDirClient=sspath
id = ROOT.EdbID(brick,0,0,0)
ss = sproc.ReadScanSet(id)
ss.Brick().SetID(brick)
zmin = ss.GetPlate(1).Z()
zmax = 0

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

vtx_cut = f'n>2&&flag==0&&vz>{zmin}&&vz<{zmax}'
# vtx_cut = f'n>2&&(flag==0||flag==3)&&vz>{zmin}&&vz<{zmax}&&abs(vx-{xcenter})<{cellsize/2}&&abs(vy-{ycenter})<{cellsize/2}'
dproc.ReadVertexTree(vertexrec, vtx_file, vtx_cut)
vertices = gAli.eVTX

#save vertices in root file
outputFile = ROOT.TFile(out_dir+"/vertex_selection_{}_{}.root".format(xbin*10, ybin*10),"RECREATE")	
outputTree = ROOT.TTree("vtx","Tree of vertices")

_brick = array('i', [0])
_cell = array('i', [0])
_cellx = array('i', [0])
_celly = array('i', [0])
_vID = array('i', [0])
_voi = array('i', [0])
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
_mom = array('f', [0])
_mom_long = array('f', [0])
_mom_max = array('f', [0])
_mom_cat = array('f', [0])
_mom_t = array('f', N*[0])

outputTree.Branch("brick", _brick, "brick/I")
outputTree.Branch("cell", _cell, "cell/I")
outputTree.Branch("cellx", _cellx, "cellx/I")
outputTree.Branch("celly", _celly, "celly/I")
outputTree.Branch("vID", _vID, "vID/I")
outputTree.Branch("voi", _voi, "voi/I")
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
outputTree.Branch("plate", _plate, "plate/I")
outputTree.Branch("nfirst", _nfirst, "nfirst[ntrks]/I")
outputTree.Branch("nlast", _nlast, "nlast[ntrks]/I")
outputTree.Branch("tTX", _tTX, "tTX[ntrks]/F")
outputTree.Branch("tTY", _tTY, "tTY[ntrks]/F")
outputTree.Branch("tLastX", _tLastX, "tLastX[ntrks]/F")
outputTree.Branch("tLastY", _tLastY, "tLastY[ntrks]/F")
outputTree.Branch("tLastZ", _tLastZ, "tLastZ[ntrks]/F")
outputTree.Branch("tID", _tID, "tID[ntrks]/I")
outputTree.Branch("prob", _prob, "prob/F")
outputTree.Branch("maxaperture", _maxaperture, "maxaperture/F")
outputTree.Branch("meanaperture", _meanaperture, "meanaperture/F")
outputTree.Branch("dphi", _dphi, "dphi[ntrks]/F")
outputTree.Branch("maxdphi", _maxdphi, "maxdphi/F")
outputTree.Branch("magdphi", _magdphi, "magdphi/F")
outputTree.Branch("meanphi", _meanphi, "meanphi/F")
outputTree.Branch("mom", _mom, "mom/F")
outputTree.Branch("mom_long", _mom_long, "mom_long/F")
outputTree.Branch("mom_max", _mom_max, "mom_max/F")
outputTree.Branch("mom_cat", _mom_cat, "mom_cat/F")
outputTree.Branch("mom_t", _mom_t, "mom_t[ntrks]/F")

vtx_of_interest = 0
for vtx in vertices:
    vtx_mom = 0
    mom_max = 0
    mom_long = 0
    seg_max = 0
    valid_rec = 0
    cat = 0

    vz = vtx.VZ()
    vx = vtx.VX()
    vy = vtx.VY()
    flag = vtx.Flag()
    ntrks = vtx.N()
    # if vtx.ID() in voiid:
    #     print(vtx.ID, ntrks)
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
        npl = track.Npl()
        tx = track.TX()
        ty = track.TY()
        sx = track.GetSegmentLast().X()
        sy = track.GetSegmentLast().Y()
        sz = track.GetSegmentLast().Z()
        nfirst = track.GetSegmentFirst().Plate()
        nlast = track.GetSegmentLast().Plate()
        plateList.append(nfirst)
        nava = from_plate - nfirst + 1
        FF = float(nseg)/float(nava)
        ffList.append(FF)
        impact_parameter = vtx.GetVTa(itrack).Imp()
        ipList.append(impact_parameter)
        TNorm = ROOT.TMath.Sqrt(tx*tx + ty*ty)
        TXList.append(track.TX())
        TYList.append(track.TY())
        TXNList.append(track.TX()/TNorm)
        TYNList.append(track.TY()/TNorm)
        phiList.append(track.Phi())
        _nseg[itrack] = nseg
        _npl[itrack] = npl
        _fillfact_t[itrack] = FF
        _ip[itrack] = impact_parameter
        _tTX[itrack] = tx
        _tTY[itrack] = ty
        _tID[itrack] = track.Track()
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
            tx= tx - t2.TX()
            ty= ty - t2.TY()
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

    _brick[0] = brick
    _cell[0] = icell
    _cellx[0] = xbin
    _celly[0] = ybin
    _vx[0]=vx
    _vy[0]=vy
    _vz[0]=vz
    _vID[0] = vtx.ID()
    # if vtx.ID() in voiid: _voi[0] = 1
    # else: _voi[0] = 0
    _ntrks[0] = ntrks
    _fillfact[0] = np.mean(ffList)
    _meanIP[0] = np.mean(ipList)
    _prob[0] = vtx.V().prob()
    _maxaperture[0] = vtx.MaxAperture()
    _maxdphi[0] = np.max(dPhiList)
    _magdphi[0] = magdphi
    _meanphi[0] = np.mean(phiList)
    _meanaperture[0] = np.mean(apeList)
    _mom[0] = vtx_mom
    _mom_long[0] = mom_long
    _mom_max[0] = mom_max
    _mom_cat[0] = cat
    outputTree.Fill()

del gAli

#write output files
outputFile.cd()
outputTree.Write()
outputFile.Close()

end_time = time()-start_time
print("TOTAL ELAPSED TIME ", end_time)
