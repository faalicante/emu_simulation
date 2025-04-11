#include <sys/stat.h>

int file_exists(const char *filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
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

void merge_events(int brickID) {
    int nentr = 4268;
    EdbID        idset;
    TObjArray    trk_out, vtx_out;
    EdbScanProc  gSproc;
    EdbVertexRec gEVR;
    gSproc.eProcDirClient=".";
    EdbScanCond  *scancond = new EdbScanCond();
    scancond->SetDegrad(3);
    scancond->SetSigma0(1.5,1.5,0.0015,0.0015);
    EdbPVRec     vtx_in;
    EdbPVRec *trk_in = new EdbPVRec();
    trk_in->SetScanCond(scancond);
    EdbDataProc *dproc = new EdbDataProc();
    vtx_in.SetScanCond(scancond);
    gEVR.eEdbTracks = vtx_in.eTracks;
    gEVR.eVTX       = vtx_in.eVTX;
    gEVR.SetPVRec(&vtx_in);
    gEVR.eDZmax=3000.;
    gEVR.eProbMin=0.001;
    gEVR.eImpMax=3.5;
    gEVR.eUseMom=0;
    gEVR.eUseSegPar=0;
    gEVR.eQualityMode=0;
    for (int evID=0; evID<nentr; evID++) {
        idset.Set(Form("%i.0.0.%i", brickID, evID+1));
        TString name;
        gSproc.MakeFileName(name,idset,"vtx.root",false);
        if (!file_exists(name.Data())) {
            std::cout << "File " << name.Data() << " does not exist. Skipping event." << std::endl;
            continue;
        }
        dproc->ReadVertexTree(gEVR, name.Data(), "1");
        gSproc.ReadTracksTree(idset, *trk_in, "1");
    }
    add_to_collection(*trk_in, vtx_in, trk_out, vtx_out);

    EdbDataProc::MakeVertexTree(vtx_out,Form("b0000%i/b0000%i.0.0.0.vtx.root", brickID, brickID));
    EdbDataProc::MakeTracksTree(trk_out,0,0,Form("b0000%i/b0000%i.0.0.0.trk.root", brickID, brickID));
}   