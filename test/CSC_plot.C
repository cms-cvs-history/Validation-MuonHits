{
std::map<int, TH1*> mh_nCSCHits;
std::map<int, TH1*> mh_nCSCMuonHits;
std::map<int, TH1*> mh_CSC_enloss;
std::map<int, TH1*> mh_CSC_tof;

std::ostringstream ss;

Int_t id, nevents=0;
Float_t pow6=1000000.0;

mh_nCSCHits.clear();
mh_nCSCMuonHits.clear();
mh_CSC_enloss.clear();
mh_CSC_tof.clear();

/// Name input data root file, output histogram root file and tree

  char * inprootfilename = "muonsimvalid.root";
  char * outrootfilename = "muonsimvalid_hist.root";
  char * treename = "Events"; 
  
  hist_file=new TFile(outrootfilename,"RECREATE");
  hist_file->cd();

  delete gROOT->GetListOfFiles()->FindObject(inprootfilename);
  TFile * myf  = new TFile(inprootfilename);  
  TTree * tree = dynamic_cast<TTree*>(myf->Get("Events"));
  assert(tree != 0);

  nevents = tree->GetEntries();
  std::cout << "Number of events = " << nevents << std::endl;

/// Choose the object to work with

  TBranch *vpbrnch = tree->GetBranch("PMuonSimHit_vp_Hits_MuonHits.obj");
  assert(vpbrnch != 0);

  PMuonSimHit vp;
  vpbrnch->SetAddress(&vp);

/// Fill the histograms

hist_file->cd();
for (Int_t ev=1; ev<=nevents; ev++) {
   vpbrnch->GetEntry( ev );

   /// Select the CSC subdetector
 
   std::vector<PMuonSimHit::CSC> CSC   = vp.getCSCHits();

   /// Number of all CSC hits

   Int_t nCSCHits   = vp.getnCSCHits();
   id=1;
   if (mh_nCSCHits.count(id) == 0) {
      ss<<"Number of all CSC hits_"<<id;
      mh_nCSCHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 100.0);
      mh_nCSCHits[id]->SetXTitle("Number of all hits in CSC");
      mh_nCSCHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nCSCHits[id]->Fill((float)nCSCHits,1.0);

   /// Number of muon hits in CSC

   Int_t nCSCMuonHits=0;
   for (Int_t i = 0; i < CSC.size(); ++i) 
     if(CSC[i]._particleType==13)  nCSCMuonHits++;   
   id=2;
   if (mh_nCSCMuonHits.count(id) == 0) {
      ss<<"Number of muon CSC hits_"<<id;
      mh_nCSCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 50.0);
      mh_nCSCMuonHits[id]->SetXTitle("Number of muon hits in CSC");
      mh_nCSCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nCSCMuonHits[id]->Fill((float)nCSCMuonHits,1.0);     

   for (Int_t i = 0; i < CSC.size(); ++i) {

   /// Select CSC muon hits only for histograms below

     if(CSC[i]._particleType==13) {

   /// Plot CSC chambers identified by endcap,station and ring only

       id=CSC[i]._cscId/1000;

   /// Energy losses in CSC

       Float_t eloss=CSC[i]._enloss*pow6;
       if (mh_CSC_enloss.count(id) == 0) {
          ss<<"ME"<<id<<"_energy_loss";
          mh_CSC_enloss[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 50.0);
          mh_CSC_enloss[id]->SetXTitle("Energy Loss(keV)");
          mh_CSC_enloss[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_CSC_enloss[id]->Fill(eloss,1.0);

   /// Time of flight for CSC

       Float_t tof=CSC[i]._tof;
       if (mh_CSC_tof.count(id) == 0) {
          ss<<"ME"<<id<<"_tof";
          mh_CSC_tof[id]  = new TH1F(ss.str().c_str(),"", 60, 0.0, 60.0);
          mh_CSC_tof[id]->SetXTitle("Time of Flight (ns)");
          mh_CSC_tof[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_CSC_tof[id]->Fill(tof,1.0);
     }
   }
}
hist_file->Write();
delete hist_file;
}
