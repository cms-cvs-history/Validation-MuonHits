{
std::map<int, TH1*> mh_nDTHits;
std::map<int, TH1*> mh_nDTMuonHits;
std::map<int, TH1*> mh_wheel;
std::map<int, TH1*> mh_station;
std::map<int, TH1*> mh_sector;
std::map<int, TH1*> mh_superlayer;
std::map<int, TH1*> mh_layer;
std::map<int, TH1*> mh_wire;
std::map<int, TH1*> mh_chambocc;
std::map<int, TH1*> mh_DT_enloss;
std::map<int, TH1*> mh_DT_mom_stat1;
std::map<int, TH1*> mh_DT_mom_stat4;
std::map<int, TH1*> mh_DT_mom_stat14;
std::map<int, TH1*> mh_DT_radius;
std::map<int, TH1*> mh_DT_costeta;
std::map<int, TH2*> mh_DT_pos;

std::ostringstream ss;

unsigned int iden;
Int_t id,  nevents=0;
Int_t wheel, station, sector, superlayer, layer, wire=100;
Int_t path, pathchamber=0;
Int_t touch1, touch4;
Float_t pow6=1000000.0;
Float_t mom1, mom4 =0.;
Float_t costeta, radius,sinteta = 0.;
Float_t xposglob, yposglob = 9999.;


mh_nDTHits.clear();
mh_nDTMuonHits.clear();
mh_wheel.clear();
mh_station.clear();
mh_sector.clear();
mh_superlayer.clear();
mh_layer.clear();
mh_wire.clear();
mh_chambocc.clear();
mh_DT_enloss.clear();
mh_DT_mom_stat1.clear();
mh_DT_mom_stat4.clear();
mh_DT_mom_stat14.clear();
mh_DT_radius.clear();
mh_DT_costeta.clear();
mh_DT_pos.clear();

/// Name input data root file, output histogram root file and tree

  char * inprootfilename = "muonsimvalid.root";
  char * outrootfilename = "muonsimvalid_DT_hist.root";
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

   /// Select the DT subdetector
 
   std::vector<PMuonSimHit::DT> DT   = vp.getDTHits();

   touch1 = 0;
   touch4 = 0;

   /// Number of all DT hits

   Int_t nDTHits   = vp.getnDTHits();
   id=1;
   if (mh_nDTHits.count(id) == 0) {
      ss << "Number of all DT hits";
      mh_nDTHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 100.0);
      mh_nDTHits[id]->SetXTitle("Number of all hits in DT");
      mh_nDTHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nDTHits[id]->Fill((float)nDTHits,1.0);

//   cout << " nDTHits " << nDTHits << " DT size " << DT.size()  << std::endl;


   /// Number of muon hits in DT

   Int_t nDTMuonHits=0;
   for (Int_t i = 0; i < DT.size(); ++i)  
     if(DT[i]._particleType==13)  nDTMuonHits++;
   id=2;
   if (mh_nDTMuonHits.count(id) == 0) {
      ss << "Number of muon DT hits";
      mh_nDTMuonHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 50.0);
      mh_nDTMuonHits[id]->SetXTitle("Number of muon hits in DT");
      mh_nDTMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nDTMuonHits[id]->Fill((float)nDTMuonHits,1.0);

   for (Int_t i = 0; i < DT.size(); ++i) {

   /// Select DT muon hits only for histograms below

     if(DT[i]._particleType==13) {

   /// Plot DT chambers identified by wheel, station, sector, superlayer,
   /// layer and wire

       iden=DT[i]._detUnitId;

       wheel = ((iden>>15) & 0x7 ) -3  ;
       station = ((iden>>22) & 0x7 ) ;
       sector = ((iden>>18) & 0xf ) ;
       superlayer = ((iden>>13) & 0x3 ) ;
       layer = ((iden>>10) & 0x7 ) ;
       wire = ((iden>>3) & 0x7f ) ;

   /// Wheel occupancy

   id=3;
   if (mh_wheel.count(id) == 0) {
      ss<<" Wheel occupancy ";
      mh_wheel[id]  = new TH1F(ss.str().c_str(),"", 10, -5.0, 5.0);
      mh_wheel[id]->SetXTitle("Wheel number");
      mh_wheel[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_wheel[id]->Fill((float)wheel,1.0);

   /// Station occupancy

   id=4;
   if (mh_station.count(id) == 0) {
      ss<<" Station occupancy ";
      mh_station[id]  = new TH1F(ss.str().c_str(),"", 6, 0., 6.0);
      mh_station[id]->SetXTitle("Station number");
      mh_station[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_station[id]->Fill((float)station,1.0);

   /// Sector occupancy

   id=5;
   if (mh_sector.count(id) == 0) {
      ss<<" Sector occupancy ";
      mh_sector[id]  = new TH1F(ss.str().c_str(),"", 20, 0., 20.);
      mh_sector[id]->SetXTitle("Sector number");
      mh_sector[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_sector[id]->Fill((float)sector,1.0);

   /// Superlayer occupancy

   id=6;
   if (mh_superlayer.count(id) == 0) {
      ss<<" SuperLayer occupancy ";
      mh_superlayer[id]  = new TH1F(ss.str().c_str(),"", 5, 0., 5.);
      mh_superlayer[id]->SetXTitle("SuperLayer number");
      mh_superlayer[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_superlayer[id]->Fill((float)superlayer,1.0);

   /// Layer occupancy

   id=7;
   if (mh_layer.count(id) == 0) {
      ss<<" Layer occupancy ";
      mh_layer[id]  = new TH1F(ss.str().c_str(),"", 6, 0., 6.);
      mh_layer[id]->SetXTitle("Layer number");
      mh_layer[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_layer[id]->Fill((float)layer,1.0);

   /// Wire occupancy

   id=8;
   if (mh_wire.count(id) == 0) {
      ss<<" Wire occupancy ";
      mh_wire[id]  = new TH1F(ss.str().c_str(),"", 100, 0., 100.);
      mh_wire[id]->SetXTitle("Wire number");
      mh_wire[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_wire[id]->Fill((float)wire,1.0);

 // Cell and layer occupancy

// MB1

   id=22;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb1_lay1_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 52, 0., 52.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 1 && superlayer == 1 && layer == 1 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=23;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb1_lay2_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 52, 0., 52.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 1 && superlayer == 1 && layer == 2 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=24;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb1_lay3_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 52, 0., 52.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 1 && superlayer == 1 && layer == 3 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=25;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb1_lay4_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 52, 0., 52.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 1 && superlayer == 1 && layer == 4 )
   mh_chambocc[id]->Fill((float)wire,1.0);

// MB2

   id=26;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb2_lay1_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 62, 0., 62.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 2 && superlayer == 1 && layer == 1 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=27;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb2_lay2_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 62, 0., 62.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 2 && superlayer == 1 && layer == 2 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=28;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb2_lay3_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 62, 0., 62.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 2 && superlayer == 1 && layer == 3 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=29;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb2_lay4_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 62, 0., 62.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 2 && superlayer == 1 && layer == 4 )
   mh_chambocc[id]->Fill((float)wire,1.0);


// MB3

   id=30;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb3_lay1_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 74, 0., 74.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 3 && superlayer == 1 && layer == 1 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=31;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb3_lay2_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 74, 0., 74.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 3 && superlayer == 1 && layer == 2 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=32;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb3_lay3_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 74, 0., 74.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 3 && superlayer == 1 && layer == 3 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=33;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb3_lay4_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 74, 0., 74.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 3 && superlayer == 1 && layer == 4 )
   mh_chambocc[id]->Fill((float)wire,1.0);


// MB4

   id=34;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb4_lay1_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 98, 0., 98.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 4 && (sector == 1 || sector == 2 || sector == 3 || 
   sector == 5 || sector == 6 || sector == 7 ) && layer == 1 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=35;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb4_lay2_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 98, 0., 98.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 4 && (sector == 1 || sector == 2 || sector == 3 || 
   sector == 5 || sector == 6 || sector == 7 ) && layer == 2 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=36;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb4_lay3_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 98, 0., 98.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 4 && (sector == 1 || sector == 2 || sector == 3 || 
   sector == 5 || sector == 6 || sector == 7 ) && layer == 3 )
   mh_chambocc[id]->Fill((float)wire,1.0);

   id=37;
   if (mh_chambocc.count(id) == 0) {
      ss<<" mb4_lay4_occ ";
      mh_chambocc[id]  = new TH1F(ss.str().c_str(),"", 98, 0., 98.);
      mh_chambocc[id]->SetXTitle("Wire number");
      mh_chambocc[id]->SetYTitle("Hits");
      ss.str("");       
   }
   if (station == 4 && (sector == 1 || sector == 2 || sector == 3 || 
   sector == 5 || sector == 6 || sector == 7 ) && layer == 4 )
   mh_chambocc[id]->Fill((float)wire,1.0);


   /// Energy losses in DT

   id = 9;
   Float_t eloss=DT[i]._enloss*pow6;
   if (mh_DT_enloss.count(id) == 0) {
      ss<<"DT energy_loss";
      mh_DT_enloss[id]  = new TH1F(ss.str().c_str(),"", 100, 0.0, 10.0);
      mh_DT_enloss[id]->SetXTitle("Energy Loss(keV)");
      mh_DT_enloss[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_DT_enloss[id]->Fill(eloss,1.0);

// Define a quantity to take into account station, splayer and layer being hit.
   path = (station-1) * 40. + superlayer * 10. + layer;

       id = 10;
       if (mh_DT_enloss.count(id) == 0) {
          ss<<" path followed by muon";
          mh_DT_enloss[id]  = new TH1F(ss.str().c_str(),"", 160, 0., 160.);
          mh_DT_enloss[id]->SetXTitle("Path followed");
          mh_DT_enloss[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_DT_enloss[id]->Fill((float)path,1.0);

// Define a quantity to take into chamber being hit.
   pathchamber = (wheel+2) * 50. + (station-1.) * 12. + sector;

       id = 11;
       if (mh_DT_enloss.count(id) == 0) {
          ss<<" chamber occupancy";
          mh_DT_enloss[id]  = new TH1F(ss.str().c_str(),"", 251, 0., 251.);
          mh_DT_enloss[id]->SetXTitle("Chamber Occupancy");
          mh_DT_enloss[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_DT_enloss[id]->Fill((float)pathchamber,1.0);

   /// Muon Momentum at MB1

       id=12;

    if (mh_DT_mom_stat1.count(id) == 0) {
          ss<<"Momentum at MB1";
          mh_DT_mom_stat1[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_DT_mom_stat1[id]->SetXTitle(" Momentum (GeV/c)");
          mh_DT_mom_stat1[id]->SetYTitle("Entries");
          ss.str("");       
    }
       if (station == 1 )
        {
         if (touch1 == 0)
         {
          mom1=DT[i]._pabs;
          touch1 = 1;
          mh_DT_mom_stat1[id]->Fill(mom1,1.0);
         }
        }

   /// Muon Momentum at MB4

       id=13;

    if (mh_DT_mom_stat4.count(id) == 0) {
          ss<<"Momentum at MB4";
          mh_DT_mom_stat4[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_DT_mom_stat4[id]->SetXTitle(" Momentum (GeV/c)");
          mh_DT_mom_stat4[id]->SetYTitle("Entries");
          ss.str("");       
    }

   /// Loss of Muon Momentum in Iron (between MB1 and MB4)
       id=14;

    if (mh_DT_mom_stat14.count(id) == 0) {
          ss<<"Loss of muon Momentum in Iron";
          mh_DT_mom_stat14[id]  = new TH1F(ss.str().c_str(),"", 80, 0.0, 40.0);
          mh_DT_mom_stat14[id]->SetXTitle(" Momentum (GeV/c)");
          mh_DT_mom_stat14[id]->SetYTitle("Entries");
          ss.str("");       
    }
       if (station == 4 )
        {
         if ( touch4 == 0)
         {
          mom4=DT[i]._pabs;
          touch4 = 1;
          mh_DT_mom_stat4[13]->Fill(mom4,1.0);
          if (touch1 == 1 )
          {
            mh_DT_mom_stat14[14]->Fill(mom1-mom4,1.0);
          }
         }
        }

   /// X-Local Coordinate vs Z-Local Coordinate

       id=15;
    if (mh_DT_pos.count(id) == 0) {
          ss<<"Local x-coord. vs local z-coord of muon hit";
          mh_DT_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -150., 150., 100, -0.8, 0.8 );
          mh_DT_pos[id]->SetXTitle(" X-coord (cm)");
          mh_DT_pos[id]->SetYTitle("Z-coord (cm)");
          ss.str("");       
    }
        mh_DT_pos[id]->Fill(DT[i]._locposx, DT[i]._locposz,1.0);

   /// X-Local Coordinate vs Y-Local Coordinate

       id=16;
    if (mh_DT_pos.count(id) == 0) {
          ss<<"local x-coord. vs local y-coord of muon hit";
          mh_DT_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -150., 150., 100, -150., 150. );
          mh_DT_pos[id]->SetXTitle(" X-local coord (cm)");
          mh_DT_pos[id]->SetYTitle("Y-local coord (cm)");
          ss.str("");       
    }
        mh_DT_pos[id]->Fill(DT[i]._locposx, DT[i]._locposy,1.0);

   /// Global Coordinates

    radius = DT[i]._globposz* ( 1.+ exp(-2.*DT[i]._globposeta) )
        / ( 1. - exp(-2.*DT[i]._globposeta) ) ;

    costeta =( 1. - exp(-2.*DT[i]._globposeta) )/( 1. + exp(-2.*DT[i]._globposeta) );

    sinteta = 2. * exp(-DT[i]._globposeta) /( 1. + exp(-2.*DT[i]._globposeta) );
    xposglob = radius*sinteta*cos(DT[i]._globposphi );
    yposglob = radius*sinteta*sin(DT[i]._globposphi );

   /// Radius of hit

       id = 19;
       if (mh_DT_radius.count(id) == 0) {
          ss<<" radius of hit ";
          mh_DT_radius[id]  = new TH1F(ss.str().c_str(),"", 100, 0., 1200.);
          mh_DT_radius[id]->SetXTitle("Radius (cm)");
          mh_DT_radius[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_DT_radius[id]->Fill(radius,1.0);

   /// Costheta of hit

       id = 20;
       if (mh_DT_costeta.count(id) == 0) {
          ss<<" costheta of hit ";
          mh_DT_costeta[id]  = new TH1F(ss.str().c_str(),"", 100, -1., 1.);
          mh_DT_costeta[id]->SetXTitle(" cos(theta) ");
          mh_DT_costeta[id]->SetYTitle("Entries");
          ss.str("");       
       }
       mh_DT_costeta[id]->Fill(costeta,1.0);

   /// Path followed by muons versus radius of hit

       id=21;
    if (mh_DT_pos.count(id) == 0) {
          ss<<" path followed by muons vs radius ";
          mh_DT_pos[id]  = new TH2F(ss.str().c_str(),"", 160, 0., 160., 100, 200., 1100. );
          mh_DT_pos[id]->SetXTitle(" Radius (cm)");
          mh_DT_pos[id]->SetYTitle("Path followed by muon ");
          ss.str("");       
    }
        mh_DT_pos[id]->Fill((float)path, radius,1.0);


   /// Z-Global Coordinate vs X-Global Coordinate

       id=17;
    if (mh_DT_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global z-coord of muon hit";
          mh_DT_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -800., 800., 100, -800., 800. );
          mh_DT_pos[id]->SetXTitle(" Global X-coord (cm)");
          mh_DT_pos[id]->SetYTitle("Global Z-coord (cm)");
          ss.str("");       
    }
        mh_DT_pos[id]->Fill( DT[i]._globposz,xposglob, 1.0);

   /// X-Global Coordinate vs Y-Global Coordinate

       id=18;
    if (mh_DT_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global y-coord of muon hit";
          mh_DT_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -800., 800., 100, -800., 800. );
          mh_DT_pos[id]->SetXTitle(" Global X-coord (cm)");
          mh_DT_pos[id]->SetYTitle("Global Y-coord (cm)");
          ss.str("");       
    }
        mh_DT_pos[id]->Fill(xposglob, yposglob,1.0);



   }
  }
}
hist_file->Write();
delete hist_file;
}
