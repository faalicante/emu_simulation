int neutrino_inbkg(int event, int cell){

  const int nplates = 60;
  int brickID = 21;
  TString cutstring = Form("s.eMCEvt==%i", event);
  TCut *cut = new TCut(cutstring);
  TString newpath = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_muon1.3E5";
  TString nupath = "/eos/experiment/sndlhc/MonteCarlo/FEDRA/numucc/numucc_eff9_smear1_evt+1";
  TString mupath = Form("/eos/experiment/sndlhc/MonteCarlo/FEDRA/muon1.3E5/cell_reco/%i", cell);
  
  for (int i = 1; i <= 60; i++) {
    EdbCouplesTree *ect = new EdbCouplesTree();
    EdbCouplesTree *ect_nu = new EdbCouplesTree();
    EdbCouplesTree *ect_mu = new EdbCouplesTree();
    TString cpname = Form("/b%06i/p%03i/%i.%i.0.0.cp.root", brickID, i, brickID, i);
    ect->InitCouplesTree("couples", newpath+cpname, "RECREATE");
    ect_nu->InitCouplesTree("couples", nupath+cpname, "READ");
    ect_mu->InitCouplesTree("couples", mupath+cpname, "READ");
    
    if (!(ect_nu->eTree)||!(ect_mu->eTree)){
      cout << "File not existing" << endl;
      return -1;
    }
    ect_nu->eCut = *cut;
    TEventList *cutlist = ect_nu->InitCutList();
    if (!cutlist){ 
      cout << "We have no entries" << endl;
      return -1;
    }
    
    // adding neutrino couples
    int nsegcut = cutlist->GetN();
    cut->Print();
    cout << "We have "<< nsegcut<< " good couples in plate "<< i << endl;
    int ihit = 0;
    for (int ientry=0; ientry < nsegcut; ientry++){
      cout << "Neutrino entry " << ientry << " " << ihit << endl;
      int iseg = cutlist->GetEntry(ientry);
      ect_nu->GetEntry(iseg);
      EdbSegP *seg = ect_nu->eS;
      ect->eS->Set(ihit,seg->X(),seg->Y(),seg->TX(),seg->TY(),1,1);
      ect->eS->SetMC(event, seg->MCTrack());
      ect->eS->SetAid(seg->Aid(0), 0);
      ect->eS->SetP(seg->P());
      ect->eS->SetVid(seg->Vid(0),0);
      ect->eS->SetW(70);
      ect->eS->SetDZem(seg->DZem());
      ect->Fill();
      ihit++;
    }
    // adding muon couples
    for (int jentry=0; jentry < ect_mu->eTree->GetEntries(); jentry++){
      if (jentry%100000==0) cout << "Muon entry " << jentry << endl;
      ect_mu->GetEntry(jentry);
      EdbSegP *seg = ect_mu->eS;
      ect->eS->Set(ihit,seg->X(),seg->Y(),seg->TX(),seg->TY(),1,0);
      ect->eS->SetMC(event, seg->MCTrack());
      ect->eS->SetAid(seg->Aid(0), 0);
      ect->eS->SetP(seg->P());
      ect->eS->SetVid(seg->Vid(0),0);
      ect->eS->SetW(70);
      ect->eS->SetDZem(seg->DZem());
      ect->Fill();
      ihit++;
    }
    cout << "We have saved " << ihit << " hits in plate " << i << endl;
    ect->Close();
  }
  return 0;
}
