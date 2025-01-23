import ROOT
import os,sys
from array import array
import rootUtils as ut
import time
from datetime import date
import random
today = date.today().strftime('%d%m%y')

def GetVetoBar(detID):
    plane = int((detID/1000)%10)
    bar = int((detID%10000)%1000)
    return plane, bar

def getOriginAndDims(nodepath):
    nav = ROOT.gGeoManager.GetCurrentNavigator()
    nav.cd(nodepath)
    O = {'X':0, 'Y':0, 'Z':0}
    D = {'X':0, 'Y':0, 'Z':0}
    N = nav.GetCurrentNode()
    S = N.GetVolume().GetShape()
    D['X'], D['Y'], D['Z'] = S.GetDX(),S.GetDY(),S.GetDZ()
    O['X'], O['Y'], O['Z'] = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
    OriginArray = array('d', O.values())
    OriginTrans = array('d', [0, 0, 0])
    nav.LocalToMaster(OriginArray, OriginTrans)
    O['X'], O['Y'], O['Z'] = OriginTrans[0],OriginTrans[1],OriginTrans[2]
    return O, D

import SndlhcGeo
#pathtoSims = "/eos/user/d/dannc/muonDis_sim/ecut1.0_z-7_2.5m_Ioni_lateFLUKA/"
pathtoSims = "/eos/experiment/sndlhc/users/dancc/MuonDIS/ecut1.0_z-7_2.5m_Ioni_latelateFLUKA/"
geofile = pathtoSims+'/muonDis_201/1/geofile_full.muonDIS-TGeant4-muonDis_201.root'
geo = SndlhcGeo.GeoInterface(geofile)
TargetNode = '/cave_1/Detector_0/volTarget_1'
TargOr, TargDim = getOriginAndDims(TargetNode)
TZ_start, TZ_end = TargOr['Z']-TargDim['Z'], TargOr['Z']+TargDim['Z']
TY_start, TY_end = TargOr['Y']-TargDim['Y'], TargOr['Y']+TargDim['Y']
TX_start, TX_end = TargOr['X']-TargDim['X'], TargOr['X']+TargDim['X']

def vtxInTarget(vtx):
    volName = ''
    try:
        node = geo.sGeo.FindNode(vtx.X(), vtx.Y(), vtx.Z()).GetMotherVolume()
        volName = node.GetName()
    except:
        return False
    if volName == 'volTarget' or volName[:5] == 'Brick' or volName[:5] == 'Scifi' or volName[:5] == 'Fiber':
        return True
    else:
        return False
    """if vtx.X() >= TX_start and vtx.X() <= TX_end:
        if vtx.Y() >= TY_start and vtx.Y() <= TY_end:
            if vtx.Z()>= TZ_start and vtx.Z()<= TZ_end:
                return True
    else:
        return False"""
    


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="single input file", required=False)
parser.add_argument("-n", "--nEvents", dest="nEvents", type=int, help="number of events to process", required=False, default=0)
parser.add_argument("-clusID", dest="ClusterID", required=False, default=0)
parser.add_argument("-procID", dest="ProcID", required=False, default=None, type=int)
options = parser.parse_args()


#pathtoPlots = '/afs/cern.ch/work/d/dannc/public/DISdata_analysis/'+str(today)+'/'
pathtoPlots = '/eos/user/d/dannc/DIS_neuanalysis/'+str(today)+'/'


if not os.path.exists(pathtoPlots):
    os.makedirs(pathtoPlots)

crsecfilepath = '/eos/experiment/sndlhc/MonteCarlo/Pythia6/MuonDIS/muDIScrossSec.root'
if not os.path.exists(crsecfilepath):
    crsecfilepath = 'muDIScrossSec.root'

start_time = time.time()

muonSample = []
for i in range(1, 11):
    muonSample.append('muonDis_2'+str(i).zfill(2))
    muonSample.append('muonDis_5'+str(i).zfill(2))

########

h               = {} #histo container
obj             = {} #objects container
MaxEvts         = 0
failFindNode 	= []

mu_start, mu_end = -700-10, 250+10

Ints      = {1:'', 2:'_intTarget', 4:'_noVeto'}
typelist = ['', '_firsttrack', '_maxentrack', '_randomtrack']
#Histos
ut.bookHist(h, 'MuEnergy', ';E_{#mu};dN/dE [GeV^{-1}]', 2*480, 0, 4800)
ut.bookHist(h, 'MuEnergy_noSND', 'Muons not interacting in SND;E_{#mu};dN/dE [GeV^{-1}]', 2*480, 0, 4800)
ut.bookHist(h, 'muVtx_XY', 'Muon interactions;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
ut.bookHist(h, 'muVtx_YZ', 'Muon interactions;Z [cm];Y [cm]', mu_end-mu_start, mu_start, mu_end, 300, -60, 240)
ut.bookHist(h, 'muVtx_XZ', 'Muon interactions;Z [cm];X [cm]', mu_end-mu_start, mu_start, mu_end, 400, -200, 200)
ut.bookHist(h, 'muVtx_XY_noSND', 'Muon interactions not in SND;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
ut.bookHist(h, 'muVtx_YZ_noSND', 'Muon interactions not in SND;Z [cm];Y [cm]', mu_end-mu_start, mu_start, mu_end, 300, -60, 240)
ut.bookHist(h, 'muVtx_XZ_noSND', 'Muon interactions not in SND;Z [cm];X [cm]', mu_end-mu_start, mu_start, mu_end, 400, -200, 200)
## Add the XZ view here
for Int in Ints.values():    
    for _track in typelist:
        ut.bookHist(h, 'Energy_neutrals'+Int+_track, 'All neutrals energy'+Int+_track+';E_0;dN/dE [5 GeV^{-1}]', 2*480, 0, 4800)
        ut.bookHist(h, 'TX_neutrals'+Int+_track, 'All neutrals TX'+Int+_track+';TX;dN', 300, -1.5, 1.5)
        ut.bookHist(h, 'TY_neutrals'+Int+_track, 'All neutrals TY'+Int+_track+';TY;dN', 300, -1.5, 1.5)
        ut.bookHist(h, 'XY_neutrals'+Int+_track, 'All neutrals XY'+Int+_track+';X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
        ut.bookHist(h, 'YZ_neutrals'+Int+_track, 'All neutrals YZ'+Int+_track+';Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)





obj['mainChain'] = ROOT.TChain('cbmsim')
obj['crsecfile'] = ROOT.TFile(crsecfilepath)
obj['g_13'] = obj['crsecfile'].Get('g_13').Clone('g_13')
obj['g_-13'] = obj['crsecfile'].Get('g_-13').Clone('g_-13')
sTree = obj['mainChain']

"""if options.ProcID is not None:
    fname="sndLHC.muonDIS-TGeant4-"+muonSample[options.ProcID]+"_digCPP.root"
    nfolders = 40
    if options.ProcID == 18 or options.ProcID == 19: nfolders = 27
    for ifile in range(1, nfolders):
        if ifile == 1:  source = ROOT.FairFileSource(pathtoSims+muonSample[options.ProcID]+'/'+str(ifile)+'/'+fname)
        else: source.AddFile(pathtoSims+muonSample[options.ProcID]+'/'+str(ifile)+'/'+fname)
        sTree.Add(pathtoSims+muonSample[options.ProcID]+'/'+str(ifile)+'/'+fname)
        print('Adding '+pathtoSims+muonSample[options.ProcID]+'/'+str(ifile)+'/'+fname+' to the chain')"""
FileDict = {}
nfiles = 0
for folder in range(20):
    fname="sndLHC.muonDIS-TGeant4-"+muonSample[folder]+"_digCPP.root"
    #nfolders=40
    #if folder == 18 or folder == 19: nfolders = 27
    nfolders = 101
    for ifile in range(1, nfolders):
        FileDict[nfiles]=pathtoSims+muonSample[folder]+'/'+str(ifile)+'/'+fname
        nfiles+=1

if options.ProcID is not None:
    sTree.Add(FileDict[options.ProcID])
    print('Reading file: '+FileDict[options.ProcID])

if options.nEvents != 0:
    MaxEvts = options.nEvents
else:
    MaxEvts = sTree.GetEntries()

if MaxEvts > sTree.GetEntries():
    MaxEvts = min(MaxEvts, sTree.GetEntries())


print('There will be processed '+str(MaxEvts)+' events!')

prefix = 'MC'
"""if options.ProcID is not None:
    tag = '.'+muonSample[options.ProcID]"""
if options.ProcID is not None:
    tag = '.'+str(options.ProcID)

outFilename = pathtoPlots+'/'+str(options.ClusterID)+'.'+prefix+'-DISneutrals'+tag+'.root'
obj['outFile'] = ROOT.TFile(outFilename, 'RECREATE')
#sTree.SetBranchStatus("EventHeader", 0)
#sTree.SetBranchStatus("EventHeader.", 0)
sTree.GetEntry(0)
if sTree.GetFile().Get("FileHeader"): sTree.GetFile().Get("FileHeader").Write()
if sTree.GetFile().Get("FileHeaderHeader"): sTree.GetFile().Get("FileHeaderHeader").Write()
sTree.GetFile().Get("BranchList").Write("BranchList", ROOT.TObject.kSingleKey)
sTree.GetFile().Get("TimeBasedBranchList").Write("TimeBasedBranchList", ROOT.TObject.kSingleKey)


outTree = sTree.CloneTree(0)

norm            = {}
eventWeight     = {}

norm[0] = 10.256410     # 1.8x1.8 m2 scoring plane normalization 10e-5 fb-1
norm[1] = 80.
norm[2] = 8e8/5e7
norm[3] = 8e8/2e8
nMult = 10
failFindNode 	= []

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()

Neutrals_list = [2112, -2112, 310, 130, 311, 22]
# Event loop
for i_event, event in enumerate(sTree) :
    if i_event > MaxEvts : break
    if i_event%1000 == 0: print("Sanity check, current event ", i_event)
    muontrack = event.MCTrack[0]
    Mupdg = muontrack.GetPdgCode()
    # Normalizations
    w       = norm[2]*muontrack.GetWeight()
    wLHC    = w/nMult/2     # I am using the FLUKA sample 2 times for mu->p and mu->n
    wInter  = event.MCTrack[2].GetWeight()
    wDis    = 0.6e-3*obj['g_'+str(Mupdg)].Eval(muontrack.GetEnergy())
    finalWEIGHT  = wLHC*wInter*wDis
    eventWeight[i_event] = finalWEIGHT

    writeEvt = False

    #if len(event.Digi_ScifiHits) == 0: continue     # events with Digitised ScifiHits

    intTarget       = False
    intElsewhere    = False
    intMuFilter     = False
    #volName = ''
    obj['muVtx'] = ROOT.TVector3(muontrack.GetStartX(), muontrack.GetStartY(), muontrack.GetStartZ())
    muVtx = obj['muVtx']
    """
    try:
        obj['muIntVol'] = geo.sGeo.FindNode(muVtx.X(), muVtx.Y(), muVtx.Z()).GetMotherVolume()
        volName = obj['muIntVol'].GetName()
    except:
        failFindNode.append(i_event)
        continue
    if volName == 'volTarget' or volName[:5] == 'Brick' or volName[:5] == 'Scifi' or volName[:5] == 'Fiber':
        intTarget = True
    elif volName == 'volMuFilter' or volName[:7] == 'volMuUp' or volName[:10] == 'volFeBlock' or volName[:9] == 'volMuDown' or volName[:8] == 'subDSBox' or volName[:8] == 'subUSBox' or volName == 'volVeto' or volName[:7] == 'volVeto' or volName[:10] == 'subVetoBox':
        intMuFilter = True
    else:
        intElsewhere = True
    """
    h['MuEnergy'].Fill(muontrack.GetEnergy(), eventWeight[i_event])
    h['muVtx_XY'].Fill(muVtx.X(), muVtx.Y(), eventWeight[i_event])
    h['muVtx_YZ'].Fill(muVtx.Z(), muVtx.Y(), eventWeight[i_event])
    h['muVtx_XZ'].Fill(muVtx.Z(), muVtx.X(), eventWeight[i_event])
    #if intTarget or intMuFilter: continue
    #h['MuEnergy_noSND'].Fill(muontrack.GetEnergy(), eventWeight[i_event])
    #h['muVtx_XY_noSND'].Fill(muVtx.X(), muVtx.Y(), eventWeight[i_event])
    #h['muVtx_YZ_noSND'].Fill(muVtx.Z(), muVtx.Y(), eventWeight[i_event])


    VetoFired = False
    n_neutrals = 0
    nlist = list()
    
    for i_track, track in enumerate(event.MCTrack):
        MotherID = track.GetMotherId()
        Energy = track.GetEnergy()
        PdgCode = track.GetPdgCode()
        if Energy < 5.: continue
        #if PdgCode not in Neutrals_list: continue

        #Get neutrals
        if not MyPDG.GetParticle(PdgCode):
            if PdgCode not in failedPDGs: failedPDGs.append(PdgCode)
            continue
        if MyPDG.GetParticle(PdgCode).Charge() != 0: continue
        #if ROOT.TMath.Abs(PdgCode) not in Neutrals_list: continue
        neu_vtx = ROOT.TVector3(track.GetStartX(), track.GetStartY(), track.GetStartZ())
        if vtxInTarget(neu_vtx): continue #exclude neutrals starting in the target
        TX, TY = None, None
        if track.GetPz() != 0:
            TX = track.GetPx()/track.GetPz()
            TY = track.GetPy()/track.GetPz()
        
        neutralInTarget = False
        hasScifiHit = False
        digitisedHit = False
        h['Energy_neutrals'].Fill(Energy, eventWeight[i_event])
        h['TX_neutrals'].Fill(TX, eventWeight[i_event])
        h['TY_neutrals'].Fill(TY, eventWeight[i_event])
        h['XY_neutrals'].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
        h['YZ_neutrals'].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])
        Neu_name = MyPDG.GetParticle(PdgCode).GetName()
        if 'Energy_'+Neu_name not in h.keys():
            ut.bookHist(h, 'Energy_'+Neu_name, Neu_name+' energy;E_0;dN/dE [10 GeV^{-1}]', 2*480, 0, 4800)
            ut.bookHist(h, 'TX_'+Neu_name, Neu_name+' TX;TX;dN', 300, -1.5, 1.5)
            ut.bookHist(h, 'TY_'+Neu_name, Neu_name+' TY;TY;dN', 300, -1.5, 1.5)
            ut.bookHist(h, 'XY_'+Neu_name, Neu_name+' XY;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
            ut.bookHist(h, 'YZ_'+Neu_name, Neu_name+' YZ;Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)
        h['Energy_'+Neu_name].Fill(Energy, eventWeight[i_event])
        h['TX_'+Neu_name].Fill(TX, eventWeight[i_event])
        h['TY_'+Neu_name].Fill(TY, eventWeight[i_event])
        h['XY_'+Neu_name].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
        h['YZ_'+Neu_name].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])


        # Now find interaction points
        for j_track, track2 in enumerate(event.MCTrack):
            if track2.GetMotherId()!= i_track: continue         #get daughter from neutral
            #if ROOT.TMath.Abs(track2.GetPdgCode()) in Neutrals_list: continue  
            #if not MyPDG.GetParticle(track2.GetPdgCode()): continue
            #if MyPDG.GetParticle(track2.GetPdgCode()).Charge() == 0: continue # get charged daughter
            vtx = ROOT.TVector3(track2.GetStartX(), track2.GetStartY(), track2.GetStartZ())
            if not vtxInTarget(vtx): continue
            neutralInTarget = True
            writeEvt= True
            break
                       

        if neutralInTarget:
            n_neutrals+=1
            nlist.append([i_track, Energy])
            """h['Energy_neutrals_intTarget'].Fill(Energy, eventWeight[i_event])
            h['TX_neutrals_intTarget'].Fill(TX, eventWeight[i_event])
            h['TY_neutrals_intTarget'].Fill(TY, eventWeight[i_event])
            h['XY_neutrals_intTarget'].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
            h['YZ_neutrals_intTarget'].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])"""
            for aHit in event.Digi_MuFilterHits:
                if not aHit.isValid(): continue
                if aHit.GetSystem()==1:
                    VetoFired = True
                    break
            """if not VetoFired:
                h['Energy_neutrals_noVeto'].Fill(Energy, eventWeight[i_event])
                h['TX_neutrals_noVeto'].Fill(TX, eventWeight[i_event])
                h['TY_neutrals_noVeto'].Fill(TY, eventWeight[i_event])
                h['XY_neutrals_noVeto'].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
                h['YZ_neutrals_noVeto'].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])"""


    
    if writeEvt:
        maxpair = max(nlist,key=lambda item:item[1])
        ntracks = [itrack for itrack, en in nlist]
        _MaxEnTrack = maxpair[0]
        #_rndmTrack = random.choice(ntracks)
        MaxEnTrack = event.MCTrack[_MaxEnTrack]
        FirstNeutralTrack = event.MCTrack[min(ntracks)]
        #RandomNTrack = event.MCTrack[_rndmTrack]
        """print('Ev ', i_event, ' I have ', n_neutrals, ' neutrals: ', nlist)
        print('        MaxEnTrack', _MaxEnTrack, ' FirstNeutralTrack ', min(ntracks), ' RandomNTrack ', _rndmTrack)
        print('')"""
        Tracktypes = {'_firsttrack': FirstNeutralTrack, '_maxentrack': MaxEnTrack}#, '_randomtrack': RandomNTrack}
        for _type, t in Tracktypes.items():
            pname = MyPDG.GetParticle(t.GetPdgCode()).GetName()

            if 'Energy_neutrals'+_type+'_'+pname not in h.keys():
                ut.bookHist(h, 'Energy_neutrals'+_type+'_'+pname,  pname+' energy'+_type+';E;dN/dE [10 GeV^{-1}]', 2*480, 0, 4800)
                ut.bookHist(h, 'TX_neutrals'+_type+'_'+pname, pname+'TX'+_type+';TX;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'TY_neutrals'+_type+'_'+pname, pname+'TY'+_type+';TY;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'XY_neutrals'+_type+'_'+pname, pname+'XY'+_type+';X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
                ut.bookHist(h, 'YZ_neutrals'+_type+'_'+pname, pname+'YZ'+_type+';Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)
                ut.bookHist(h, 'Energy_neutrals_noVeto'+_type+'_'+pname,  pname+' energy'+_type+' noVeto;E;dN/dE [10 GeV^{-1}]', 2*480, 0, 4800)
                ut.bookHist(h, 'TX_neutrals_noVeto'+_type+'_'+pname, pname+'TX'+_type+' noVeto;TX;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'TY_neutrals_noVeto'+_type+'_'+pname, pname+'TY'+_type+' noVeto;TY;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'XY_neutrals_noVeto'+_type+'_'+pname, pname+'XY'+_type+' noVeto;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
                ut.bookHist(h, 'YZ_neutrals_noVeto'+_type+'_'+pname, pname+'YZ'+_type+' noVeto;Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)

            h['Energy_neutrals'+_type+'_'+pname].Fill(t.GetEnergy(), eventWeight[i_event])
            h['TX_neutrals'+_type+'_'+pname].Fill(t.GetPx()/t.GetPz(), eventWeight[i_event])
            h['TY_neutrals'+_type+'_'+pname].Fill(t.GetPy()/t.GetPz(), eventWeight[i_event])
            h['XY_neutrals'+_type+'_'+pname].Fill(t.GetStartX(), t.GetStartY(), eventWeight[i_event])
            h['YZ_neutrals'+_type+'_'+pname].Fill(t.GetStartZ(), t.GetStartY(), eventWeight[i_event])

            h['Energy_neutrals'+_type].Fill(t.GetEnergy(), eventWeight[i_event])
            h['TX_neutrals'+_type].Fill(t.GetPx()/t.GetPz(), eventWeight[i_event])
            h['TY_neutrals'+_type].Fill(t.GetPy()/t.GetPz(), eventWeight[i_event])
            h['XY_neutrals'+_type].Fill(t.GetStartX(), t.GetStartY(), eventWeight[i_event])
            h['YZ_neutrals'+_type].Fill(t.GetStartZ(), t.GetStartY(), eventWeight[i_event])
            if not VetoFired:
                h['Energy_neutrals_noVeto'+_type+'_'+pname].Fill(t.GetEnergy(), eventWeight[i_event])
                h['TX_neutrals_noVeto'+_type+'_'+pname].Fill(t.GetPx()/t.GetPz(), eventWeight[i_event])
                h['TY_neutrals_noVeto'+_type+'_'+pname].Fill(t.GetPy()/t.GetPz(), eventWeight[i_event])
                h['XY_neutrals_noVeto'+_type+'_'+pname].Fill(t.GetStartX(), t.GetStartY(), eventWeight[i_event])
                h['YZ_neutrals_noVeto'+_type+'_'+pname].Fill(t.GetStartZ(), t.GetStartY(), eventWeight[i_event])

                h['Energy_neutrals_noVeto'+_type].Fill(t.GetEnergy(), eventWeight[i_event])
                h['TX_neutrals_noVeto'+_type].Fill(t.GetPx()/t.GetPz(), eventWeight[i_event])
                h['TY_neutrals_noVeto'+_type].Fill(t.GetPy()/t.GetPz(), eventWeight[i_event])
                h['XY_neutrals_noVeto'+_type].Fill(t.GetStartX(), t.GetStartY(), eventWeight[i_event])
                h['YZ_neutrals_noVeto'+_type].Fill(t.GetStartZ(), t.GetStartY(), eventWeight[i_event])
        outTree.Fill()


    ###################################
print('Analysed '+str(MaxEvts)+' events')

for hist in h.values():
    hist.Write()

obj['outFile'].Write()
obj['outFile'].Close()
print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')
