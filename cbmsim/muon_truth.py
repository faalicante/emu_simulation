# number of tracks in neutrino and muon vertices: You can loop through MCTracks and cut on the Parent ID.
# number of neutral vertices in muon interactions: My approach would be to first look for neutral particles with the muon parent, 
#       and then loop again through the MCTracks and find particles whose parent ID matches one of those neutrals.
#       Devi farti un loop sui carichi e per ognuo vedi se la parent è neutra




import ROOT
import os,sys
from array import array
import rootUtils as ut
import time
from datetime import date
import random
today = date.today().strftime('%d%m%y')
import SndlhcGeo


pathSim= "/eos/experiment/sndlhc/users/dancc/PassingMu/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_z289.374023_BRICK11/9790422/"
geofile = pathSim+"geofile_full.Ntuple-TGeant4.root"
geo = SndlhcGeo.GeoInterface(geofile)
inputFile = pathSim+"sndLHC.Ntuple-muBkg1e5cm2_B11.root"
outputFileName = 'muon_truth.root'

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

if not os.path.exists(pathPlots):
    os.makedirs(pathPlots)

start_time = time.time()


h               = {} #histo container
mu_start, mu_end = -700-10, 250+10

#Histos
ut.bookHist(h, 'NProns', 'Number of prons;n', 50, 1, 50)
ut.bookHist(h, 'Energy_neutrals', 'All neutrals energy;E_0;dN/dE [5 GeV^{-1}]', 2*480, 0, 4800)
ut.bookHist(h, 'TX_neutrals', 'All neutrals TX;TX;dN', 300, -1.5, 1.5)
ut.bookHist(h, 'TY_neutrals', 'All neutrals TY;TY;dN', 300, -1.5, 1.5)
ut.bookHist(h, 'XY_neutrals', 'All neutrals XY;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
ut.bookHist(h, 'YZ_neutrals', 'All neutrals YZ;Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)



sTree = inputFile.cbmsim
outputFile = ROOT.TFile(outputFileName, 'RECREATE')

norm            = {}
eventWeight     = {}

norm[0] = 10.256410     # 1.8x1.8 m2 scoring plane normalization 10e-5 fb-1
norm[1] = 80.
norm[2] = 8e8/5e7    #8e8 collision rate lhc = fb-1 * cross_sec 5e7 pp collision in sim
norm[3] = 8e8/2e8
nMult = 10

MyPDG = ROOT.TDatabasePDG.Instance()
failedPDGs = list()

# Event loop
for i_event, event in enumerate(sTree) :
    if i_event%1000 == 0: print("Sanity check, current event ", i_event)
    muontrack = event.MCTrack[0]
    Mupdg = muontrack.GetPdgCode()
    # Normalizations
    wLHC       = norm[3]*muontrack.GetWeight()
    eventWeight[i_event] = wLHC

    writeEvt = False
    #volName = ''
    muVtx = ROOT.TVector3(muontrack.GetStartX(), muontrack.GetStartY(), muontrack.GetStartZ())
    intTarget = vtxInTarget(muVtx)
    if not intTarget: continue
    
    h['MuEnergy'].Fill(muontrack.GetEnergy(), eventWeight[i_event])
    h['muVtx_XY'].Fill(muVtx.X(), muVtx.Y(), eventWeight[i_event])
    h['muVtx_YZ'].Fill(muVtx.Z(), muVtx.Y(), eventWeight[i_event])
    h['muVtx_XZ'].Fill(muVtx.Z(), muVtx.X(), eventWeight[i_event])

    for i_track, track in enumerate(event.MCTrack):
        MotherID = track.GetMotherId()
        Energy = track.GetEnergy()
        PdgCode = track.GetPdgCode()
        if Energy < 5.: continue

        #Get neutrals
        if not MyPDG.GetParticle(PdgCode):
            if PdgCode not in failedPDGs: failedPDGs.append(PdgCode)
            continue
        charge =  MyPDG.GetParticle(PdgCode).Charge()
        if charge != 0: continue
        neu_vtx = ROOT.TVector3(track.GetStartX(), track.GetStartY(), track.GetStartZ())
        if not vtxInTarget(neu_vtx): continue #exclude neutrals starting in the target
        TX, TY = None, None
        if track.GetPz() != 0:
            TX = track.GetPx()/track.GetPz()
            TY = track.GetPy()/track.GetPz()
        
        neutralInTarget = False


        # Now find interaction points
        nprons=0
        for j_track, track2 in enumerate(event.MCTrack):
            if track2.GetMotherId()!= i_track: continue         #get daughter from neutral
            #if ROOT.TMath.Abs(track2.GetPdgCode()) in Neutrals_list: continue  
            #if not MyPDG.GetParticle(track2.GetPdgCode()): continue
            if MyPDG.GetParticle(track2.GetPdgCode()).Charge() == 0: continue # get charged daughter
            vtx = ROOT.TVector3(track2.GetStartX(), track2.GetStartY(), track2.GetStartZ())
            if not vtxInTarget(vtx): continue
            neutralInTarget = True
            nprons+=1

        if neutralInTarget:
            h['Energy_neutrals'].Fill(Energy, eventWeight[i_event])
            h['NProns'].Fill(nprons, eventWeight[i_event])
            h['TX_neutrals'].Fill(TX, eventWeight[i_event])
            h['TY_neutrals'].Fill(TY, eventWeight[i_event])
            h['XY_neutrals'].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
            h['YZ_neutrals'].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])
            Neu_name = MyPDG.GetParticle(PdgCode).GetName()
            if 'Energy_'+Neu_name not in h.keys():
                ut.bookHist(h, 'Energy_'+Neu_name, Neu_name+' energy;E_0;dN/dE [10 GeV^{-1}]', 2*480, 0, 4800)
                ut.bookHist(h, 'NProns_'+Neu_name, Neu_name+' number of prons'+Int+';n', 50, 1, 50)
                ut.bookHist(h, 'TX_'+Neu_name, Neu_name+' TX;TX;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'TY_'+Neu_name, Neu_name+' TY;TY;dN', 300, -1.5, 1.5)
                ut.bookHist(h, 'XY_'+Neu_name, Neu_name+' XY;X [cm];Y [cm]', 400, -200, 200, 300, -60, 240)
                ut.bookHist(h, 'YZ_'+Neu_name, Neu_name+' YZ;Z [cm];Y [cm]', mu_end-mu_start+300, mu_start, mu_end+300, 300, -60, 240)
            h['Energy_'+Neu_name].Fill(Energy, eventWeight[i_event])
            h['NProns_'+Neu_name].Fill(nprons, eventWeight[i_event])
            h['TX_'+Neu_name].Fill(TX, eventWeight[i_event])
            h['TY_'+Neu_name].Fill(TY, eventWeight[i_event])
            h['XY_'+Neu_name].Fill(track.GetStartX(), track.GetStartY(), eventWeight[i_event])
            h['YZ_'+Neu_name].Fill(track.GetStartZ(), track.GetStartY(), eventWeight[i_event])


outputFile.cd()          
for hist in h.values():
    hist.Write()
outputFile.Write()
outputFile.Close()

print('Done')
print('Elapsed time: '+str((time.time()-start_time)/60.)+' mins')