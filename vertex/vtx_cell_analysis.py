import ROOT as r
import numpy as np
import fedrarootlogon
from time import time
import os.path
from array import array


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-x", dest="xcell", help="xcell", required=False, type=int)
parser.add_argument("-y", dest="ycell", help="ycell", required=False, type=int)
parser.add_argument("--cell", dest="cellID", required=True, default=None, type=int)
options = parser.parse_args()

brickID = 121
from_plate = 57
to_plate = 1
ncellsx = 18
ncellsy = 18
ncellminx = 5000
ncellmaxx= 185000
ncellminy = 5000
ncellmaxy= 185000
cellsize= (ncellmaxx-ncellminx)/ncellsx
zmin = -75600.00

path = '/eos/experiment/sndlhc/emulsionData/2022/emureco_Napoli/RUN1/b{:06d}/cells'.format(brickID)
out_dir = '/afs/cern.ch/work/f/falicant/public/RUN1/b121'

icell = options.cellID
xbin=((icell // ncellsy)) + 1
ybin=((icell % ncellsy)) + 1
# xbin=options.xcell
# ybin=options.ycell

with open(out_dir+"/vertices_of_interest.txt", "r") as file:
    lines = file.readlines()[1:]

voiid = []

for line in lines:
    values = line.split()
    if int(values[1])!=xbin or int(values[2])!=ybin: continue
    voiid.append(int(values[0]))

#histo setup
h_n = r.TH1D('n','Multiplicity;multiplicity', 30, 0, 30)
h_flag = r.TH1D('flag','Flag;flag', 6, 0, 6)
h_vz = r.TH1D('vz','Vertex z position;vz[um]', 400, -80000, 5000)
h_vxy = r.TH2D('vxy', 'Vertex xy map;vx[um];vy[um]', ncellsx*10, ncellminx, ncellmaxx, ncellsy*10, ncellminy, ncellmaxy)
h_n0 = r.TH1D('n0', 'Multiplicity;multiplicity', 27, 3, 30)
h_nseg = r.TH1D('nseg', 'Number of segments;nseg', 56, 4, 60)
h_npl = r.TH1D('npl', ' Number of crossing films;npl', 56, 4, 60)
h_ff = r.TH1D('ff', 'Fill Factor;FF', 55, 0, 1.1)
h_ip = r.TH1D('ip', 'Impact parameter;ip[um]', 100, 0, 20)
h_meanff = r.TH1D('meanff', 'Mean Fill Factor;FF', 22, 0, 1.05)
h_meanip = r.TH1D('meanip', 'Mean impact parameter;ip[um]', 100, 0, 20)
h_prob = r.TH1D('prob', 'Probability;prob', 55, 0, 1.1)
h_maxape = r.TH1D('maxape', 'Max aperture;max_ape', 100, 0, 1)
h_meanape = r.TH1D('meanape', 'Mean aperture;mean_ape', 100, 0, 1)
h_meanphi = r.TH1D('meanphi', 'Mean phi;mean_phi', 160, -4, 4)
h_maxdphi = r.TH1D('maxdphi', 'Max phi diff;max_dphi', 80, 0, 4)

N=30
start_time = time()
cell_path1 = path+'/cell_{}_{}_1x1cm/b000021'.format(xbin*10, ybin*10)
cell_path2 = path+'/cell_{}_{}_1x1cm_done/b000021'.format(xbin*10, ybin*10)
#reading vertices
if os.path.isfile(cell_path1+'/b000021.0.0.0.vtx.root'):
    vtx_file = cell_path1+'/b000021.0.0.0.vtx.root'
elif os.path.isfile(cell_path2+'/b000021.0.0.0.vtx.root'):
    vtx_file = cell_path2+'/b000021.0.0.0.vtx.root'
else:
    print("Vertex file not found, ", cell_path1)
    exit(1)

#save vertices in root file
outputFile = r.TFile(out_dir+"/vertex_selection_{}_{}.root".format(xbin*10, ybin*10),"RECREATE")	
outputTree = r.TTree("vtx","Tree of vertices")

_cellx = array('i', [0])
_celly = array('i', [0])
_vID = array('i', [0])
_voi = array('i', [0])
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
_dphi = array('f', N*[0])
_maxdphi = array('f', [0])
_meanphi = array('f', [0])
_nfirst = array('i', N*[0])
_nlast = array('i', N*[0])
_tTX = array('f', N*[0])
_tTY = array('f', N*[0])
_tLastX = array('f', N*[0])
_tLastY = array('f', N*[0])
_tLastZ = array('f', N*[0])
_tID = array('i', N*[0])

outputTree.Branch("cellx", _cellx, "cellx/I")
outputTree.Branch("celly", _celly, "celly/I")
outputTree.Branch("vID", _vID, "vID/I")
outputTree.Branch("voi", _voi, "voi/I")
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
outputTree.Branch("meanphi", _meanphi, "meanphi/F")

print("opening file: ",vtx_file)
dproc = r.EdbDataProc()
gAli = dproc.PVR()
scancond = r.EdbScanCond()
scancond.SetSigma0(1.5,1.5,0.0015,0.0015) #change sigma0
scancond.SetDegrad(3) #change angular degradation
gAli.SetScanCond(scancond)
vertexrec = r.EdbVertexRec()
vertexrec.SetPVRec(gAli)
vertexrec.eDZmax=3000.
vertexrec.eProbMin=0.001
vertexrec.eImpMax=3.5
vertexrec.eUseMom=False
vertexrec.eUseSegPar=True
vertexrec.eQualityMode=0
proc = r.EdbDataProc()
dproc.ReadVertexTree(vertexrec, vtx_file, "1")
vertices = gAli.eVTX

vtx_of_interest = 0
for vtx in vertices:
    vz = vtx.VZ()
    vx = vtx.VX()
    vy = vtx.VY()
    h_vz.Fill(vz)
    if vz < zmin: continue
    if vz > 0: continue
    h_vxy.Fill(vx, vy)
    if vx < (xbin-1)*cellsize+ncellminx: continue
    if vx >= (xbin)*cellsize+ncellminx: continue
    if vy < (ybin-1)*cellsize+ncellminy: continue
    if vy >= (ybin)*cellsize+ncellminy: continue
    flag = vtx.Flag()
    h_flag.Fill(flag)
    if flag != 0 and flag != 3: continue
    ntrks = vtx.N()
    h_n.Fill(ntrks)
    if ntrks < 3: continue
    h_n0.Fill(ntrks)
    nplList = []
    nsegList = []
    apeList = []
    phiList = []
    TXList = []
    TYList = []
    ffList = []
    ipList = []
    segidx = 0

    ### TRACKS LOOP ###
    for itrack in range(ntrks):
        track = vtx.GetTrack(itrack)
        h_ip.Fill(vtx.GetVTa(itrack).Imp())         
        ipList.append(vtx.GetVTa(itrack).Imp())
        nseg = track.N()
        npl = track.Npl()
        tx = track.TX()
        ty = track.TY()
        sx = track.GetSegmentLast().X()
        sy = track.GetSegmentLast().Y()
        sz = track.GetSegmentLast().Z()
        nfirst = track.GetSegmentFirst().Plate()
        nlast = track.GetSegmentLast().Plate()
        nava = from_plate - nfirst + 1
        FF = float(nseg)/float(nava)
        nsegList.append(nseg)
        nplList.append(npl)
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
        _ip[itrack] = vtx.GetVTa(itrack).Imp()
        _tTX[itrack] = tx
        _tTY[itrack] = ty
        _tID[itrack] = track.Track()
        _nfirst[itrack] = nfirst
        _nlast[itrack] = nlast
        _tLastX[itrack] = sx
        _tLastY[itrack] = sy           
        _tLastZ[itrack] = sz           

        for jtrack in range(itrack+1, ntrks):
            t2 = vtx.GetTrack(jtrack)
            tx= track.TX() - t2.TX()
            ty= track.TY() - t2.TY()
            apeList.append(r.TMath.Sqrt( tx*tx+ty*ty ))

    h_prob.Fill(vtx.V().prob())
    h_maxape.Fill(vtx.MaxAperture())
    h_meanff.Fill(np.mean(ffList))
    h_meanip.Fill(np.mean(ipList))
    h_meanape.Fill(np.mean(apeList))
    h_meanphi.Fill(np.mean(phiList))

    arrPhi = np.array(phiList)
    arrTX = np.array(TXList)
    arrTY = np.array(TYList)
    dPhiList = []
    for itrack in range(ntrks):
        average_of_rest = r.TMath.ATan2(np.sum(arrTY) - arrTY[itrack], np.sum(arrTX) - arrTX[itrack])
        difference = r.TMath.Abs(arrPhi[itrack] - average_of_rest)
        if difference > r.TMath.Pi():
            difference = 2*r.TMath.Pi() - difference
        dPhiList.append(difference)
        _dphi[itrack]=difference
    h_maxdphi.Fill(np.max(dPhiList))

    _cellx[0] = xbin
    _celly[0] = ybin
    _vx[0]=vx
    _vy[0]=vy
    _vz[0]=vz
    _vID[0] = vtx.ID()
    if vtx.ID() in voiid: _voi[0] = 1
    else: _voi[0] = 0
    _ntrks[0] = ntrks
    _nsegtot[0] = segidx
    _fillfact[0] = np.mean(ffList)
    _meanIP[0] = np.mean(ipList)
    _prob[0] = vtx.V().prob()
    _maxaperture[0] = vtx.MaxAperture()
    _maxdphi[0] = np.max(dPhiList)
    _meanphi[0] = np.mean(phiList)
    _meanaperture[0] = np.mean(apeList)
    outputTree.Fill()

del gAli

#write output files
outputFile.cd()
outputTree.Write()
outputFile.Close()

histoFile = r.TFile(out_dir+"/hist_out_{}_{}.root".format(xbin, ybin), "RECREATE")
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
histoFile.Write()
histoFile.Close()

end_time = time()-start_time
print("TOTAL ELAPSED TIME ", end_time)
