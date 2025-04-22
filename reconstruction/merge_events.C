#include <sys/stat.h>

int file_exists(const char *filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

void add_to_collection(EdbPVRec &vtx_in, TObjArray &vtx_out) {
    int nvtx = vtx_in.Nvtx();
    for(int iv=0; iv<nvtx; iv++) {
        EdbVertex *v = vtx_in.GetVertex(iv);
        vtx_out.Add(v);
    }
}

void add_to_collection(EdbPVRec &trk_in, EdbPVRec &vtx_in, TObjArray &trk_out, TObjArray &vtx_out) {
    int ntr = trk_in.Ntracks();
    for(int it=0; it<ntr; it++) {
        EdbTrackP *t = trk_in.GetTrack(it);
        trk_out.Add(t);
    }
    int nvtx = vtx_in.Nvtx();
    for(int iv=0; iv<nvtx; iv++) {
        EdbVertex *v = vtx_in.GetVertex(iv);
        vtx_out.Add(v);
    }
}

void MakeScanCondBT(EdbScanCond &cond)
{
  cond.SetSigma0(0.5,0.5,0.0025,0.0025);
  cond.SetDegrad(5);
  cond.SetBins(3, 3, 3, 3);
  cond.SetPulsRamp0(  12., 18. );
  cond.SetPulsRamp04( 12., 18. );
  cond.SetChi2Max( 6.5 );
  cond.SetChi2PMax( 6.5 );
  cond.SetChi2Mode( 3 );
  cond.SetRadX0( 5810. );
  cond.SetName("SND_basetrack");
}

void merge_events(int brickID) {
    int nentr = 4268;
    TObjArray    trk_out, vtx_out;
    EdbPVRec     vtx_in;
    EdbScanCond  gCond;
    EdbID        idset;
    EdbScanProc  gSproc;
    EdbVertexRec gEVR;
    gSproc.eProcDirClient=".";
    EdbScanCond  *scancond = new EdbScanCond();
    scancond->SetDegrad(5);
    scancond->SetSigma0(0.5,0.5,0.0025,0.0025);
    EdbPVRec *trk_in = new EdbPVRec();
    trk_in->SetScanCond(scancond);
    MakeScanCondBT(gCond);
    vtx_in.SetScanCond(new EdbScanCond(gCond));
    gEVR.eEdbTracks = vtx_in.eTracks;
    gEVR.eVTX       = vtx_in.eVTX;
    gEVR.SetPVRec(&vtx_in);
    gEVR.eDZmax=3000.;
    gEVR.eProbMin=0.001;
    gEVR.eImpMax=3.5;
    gEVR.eUseMom=false;
    gEVR.eUseSegPar=false;
    gEVR.eQualityMode=0;
    EdbDataProc *dproc = new EdbDataProc();
    for (int evID=1028; evID<1030; evID++) {
        idset.Set(Form("%i.0.0.%i", brickID, evID+1));
        TString name = Form("b0000%i/b0000%i.0.0.%i.vtx.root", brickID, brickID, evID+1);
        if (!file_exists(name.Data())) {
            std::cout << "File " << name.Data() << " does not exist. Skipping event." << std::endl;
            continue;
        }
        dproc->ReadVertexTree(gEVR, name.Data(), "1");
        gSproc.ReadTracksTree(idset, *trk_in, "1");
    }
    add_to_collection(*trk_in, vtx_in, trk_out, vtx_out);

    EdbDataProc::MakeVertexTree(vtx_out,Form("b0000%i/b0000%i.0.0.test.vtx.root", brickID, brickID));
    EdbDataProc::MakeTracksTree(trk_out,0,0,Form("b0000%i/b0000%i.0.0.test.trk.root", brickID, brickID));
}
