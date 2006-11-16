{
std::map<int, TH1*> mh_nRPCHits;
std::map<int, TH1*> mh_nRPCMuonHits;
std::map<int, TH1*> mh_RPC_enloss;
std::map<int, TH1*> mh_RPC_tof;
std::map<int, TH2*> mh_RPC_pos;

std::ostringstream ss;

// long int iden;
unsigned int iden;
Int_t id,  nevents=0;
Int_t region, ring, station, sector, layer, subsector, roll = 100;
Int_t path, pathchamber=0;
Int_t touch1, touch4;
Int_t touche1, touche4;
Float_t pow6=1000000.0;
Float_t mom1, mom4 =0.;
Float_t mome1, mome4 =0.;
Float_t costeta, radius,sinteta = 0.;
Float_t xposglob, yposglob = 9999.;


mh_nRPCHits.clear();
mh_nRPCMuonHits.clear();
mh_RPC_enloss.clear();
mh_RPC_tof.clear();
mh_RPC_pos.clear();

/// Name input data root file, output histogram root file and tree

  char * inprootfilename = "muonsimvalid.root";
  char * outrootfilename = "muonsimvalid_RPC_hist.root";
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

   /// Select the RPC subdetector
 
   std::vector<PMuonSimHit::RPC> RPC   = vp.getRPCHits();

   touch1 = 0;
   touch4=0;
   touche1 = 0;
   touche4 = 0;

   /// Number of all RPC hits

   Int_t nRPCHits   = vp.getnRPCHits();
   id=1;
   if (mh_nRPCHits.count(id) == 0) {
      ss << "Number of all RPC hits";
      mh_nRPCHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 100.0);
      mh_nRPCHits[id]->SetXTitle("Number of all hits in RPC");
      mh_nRPCHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nRPCHits[id]->Fill((float)nRPCHits,1.0);


   /// Number of muon hits in RPC

   Int_t nRPCMuonHits=0;
   for (Int_t i = 0; i < RPC.size(); ++i) {
     if(RPC[i]._particleType==13)  nRPCMuonHits++;
   }
   id=2;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss << "Number of muon RPC hits";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 50.0);
      mh_nRPCMuonHits[id]->SetXTitle("Number of muon hits in RPC");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nRPCMuonHits[id]->Fill((float)nRPCMuonHits,1.0);

   for (Int_t i = 0; i < RPC.size(); ++i) {

   /// Select RPC muon hits only for histograms below

     if(RPC[i]._particleType==13) {

   /// Plot RPC chambers identified by region, ring, station, sector, subsector,
   /// layer and roll

       iden=RPC[i]._detUnitId;

       region = ( ((iden>>0) & 0X3) -1 )  ;
       ring = ((iden>>2) & 0X7 ) ;

       if ( ring < 3 )
       {
        if ( region == 0 ) cout << "Region - Ring inconsistency" << endl;
        ring += 1 ;
       } else {
        ring -= 5 ;
       }

       station =  ( ((iden>>5) & 0X3) + 1 ) ;
       sector =  ( ((iden>>7) & 0XF) + 1 ) ;
       layer = ( ((iden>>11) & 0X1) + 1 ) ;
       subsector =  ( ((iden>>12) & 0X7) + 1 ) ;   //  ! Ojo que el mask figura 0x7 !!
       roll =  ( (iden>>15) & 0X7)  ;


   //    if ( ev < 20 ) {
   //      std::cout << "iden " << iden << std::endl;
   //      std::cout << "RPC[i]._detUnitId " << RPC[i]._detUnitId << std::endl;
    //     printf( "RPCmuon[i]._detUnitId  %ld ",  RPC[i]._detUnitId );
   //      std::cout << "region, ring, station,sector,lay, subsect, roll " << region << ring << 
   //   station << sector << layer << subsector << roll << std::endl;
   //    }

   /// Region occupancy

   id=3;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Region occupancy ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 6, -3.0, 3.0);
      mh_nRPCMuonHits[id]->SetXTitle("Endcap/Barrel region");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   mh_nRPCMuonHits[id]->Fill((float)region,1.0);


   /// Ring occupancy

    // Barrel
   id=4;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Ring occupancy (barrel) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 8, -3., 5.0);
      mh_nRPCMuonHits[id]->SetXTitle("Ring number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_nRPCMuonHits[id]->Fill((float)ring,1.0);

   
     // Endcap
   id=26;
   if (mh_nRPCMuonHits.count(id) == 0) {   
      ss<<" Ring occupancy (endcaps) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 8, -3., 5.0);
      mh_nRPCMuonHits[id]->SetXTitle("Ring number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_nRPCMuonHits[id]->Fill((float)ring,1.0);


   /// Station occupancy

   id=5;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Station occupancy (barrel) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 8, 0., 8.);
      mh_nRPCMuonHits[id]->SetXTitle("Station number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_nRPCMuonHits[id]->Fill((float)station,1.0);
   

   id=28;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Station occupancy (endcaps) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 8, 0., 8.);
      mh_nRPCMuonHits[id]->SetXTitle("Station number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_nRPCMuonHits[id]->Fill((float)station,1.0);
 

   /// Sector occupancy

   id=6;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Sector occupancy (barrel) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 16, 0., 16.);
      mh_nRPCMuonHits[id]->SetXTitle("Sector number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_nRPCMuonHits[id]->Fill((float)sector,1.0);

   id=27;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Sector occupancy (endcaps) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 16, 0., 16.);
      mh_nRPCMuonHits[id]->SetXTitle("Sector number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_nRPCMuonHits[id]->Fill((float)sector,1.0);

   /// Layer occupancy

   id=7;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Layer occupancy ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 4, 0., 4.);
      mh_nRPCMuonHits[id]->SetXTitle("Layer number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
    mh_nRPCMuonHits[id]->Fill((float)layer,1.0);


   /// Subsector occupancy

   id=8;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Subsector occupancy (barrel) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 10, 0., 10.);
      mh_nRPCMuonHits[id]->SetXTitle("Subsector number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_nRPCMuonHits[id]->Fill((float)subsector,1.0);

   id=35;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Subsector occupancy (endcaps) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 10, 0., 10.);
      mh_nRPCMuonHits[id]->SetXTitle("Subsector number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_nRPCMuonHits[id]->Fill((float)subsector,1.0);


   /// Roll occupancy

   id=22;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Roll occupancy (barrel) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 6, 0., 6.);
      mh_nRPCMuonHits[id]->SetXTitle("Roll number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_nRPCMuonHits[id]->Fill((float)roll,1.0);

    id=34;
   if (mh_nRPCMuonHits.count(id) == 0) {
      ss<<" Roll occupancy (endcaps ) ";
      mh_nRPCMuonHits[id]  = new TH1F(ss.str().c_str(),"", 6, 0., 6.);
      mh_nRPCMuonHits[id]->SetXTitle("Roll number");
      mh_nRPCMuonHits[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_nRPCMuonHits[id]->Fill((float)roll,1.0);

   /// Energy losses in RPC

   id = 9;
   Float_t eloss=RPC[i]._enloss*pow6;
   if (mh_RPC_enloss.count(id) == 0) {
      ss<<"RPC energy_loss (barrel) ";
      mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 10.0);
      mh_RPC_enloss[id]->SetXTitle("Energy Loss(keV)");
      mh_RPC_enloss[id]->SetYTitle("Entries");
      ss.str("");       
   }
   if (region == 0 ) mh_RPC_enloss[id]->Fill(eloss,1.0);

   id = 29;
   Float_t eloss=RPC[i]._enloss*pow6;
   if (mh_RPC_enloss.count(id) == 0) {
      ss<<"RPC energy_loss (endcaps) ";
      mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 50, 0.0, 10.0);
      mh_RPC_enloss[id]->SetXTitle("Energy Loss(keV)");
      mh_RPC_enloss[id]->SetYTitle("Entries");
      ss.str("");
   }
   if (region != 0 ) mh_RPC_enloss[id]->Fill(eloss,1.0);


// Define a quantity to take into account station, splayer and layer being hit.
 path = (region+1) * 50. + (ring+2) * 10. + (station -1) *2+ layer;
 if (region != 0) path -= 10 ; 

     id = 10;
     if (mh_RPC_enloss.count(id) == 0) {
        ss<<" path followed by muon";
        mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 160, 0., 160.);
        mh_RPC_enloss[id]->SetXTitle("Path followed");
        mh_RPC_enloss[id]->SetYTitle("Entries");
        ss.str("");
     }
     mh_RPC_enloss[id]->Fill((float)path,1.0);


   /// Muon Momentum at RB1 (Barrel)

       id=12;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Momentum at RB1";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }
  if ( region == 0 )  //  BARREL
  {
       if (station == 1 && layer == 1 )
        {
         if (touch1 == 0)
         {
          mom1=RPC[i]._pabs;
          touch1 = 1;
          mh_RPC_tof[id]->Fill(mom1,1.0);
         }
        }

   /// Muon Momentum at RB4 (Barrel)

       id=13;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Momentum at RB4";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }

   /// Loss of Muon Momentum in Iron (between RB1_layer1 and RB4)
       id=14;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Loss of muon Momentum in Iron (barrel) ";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 0.0, 40.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }
       if (station == 4 )
        {
         if ( touch4 == 0)
         {
          mom4=RPC[i]._pabs;
          touch4 = 1;
          mh_RPC_tof[13]->Fill(mom4,1.0);
          if (touch1 == 1 )
          {
            mh_RPC_tof[14]->Fill(mom1-mom4,1.0);
          }
         }
        }

  }  // End of Barrel

   /// Muon Momentum at RE1 (Endcaps)

       id=23;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Momentum at RE1";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }
  if ( region != 0 )  //  ENDCAPS
  {
       if (station == 1 )
        {
         if (touche1 == 0)
         {
          mome1=RPC[i]._pabs;
          touche1 = 1;
          mh_RPC_tof[id]->Fill(mome1,1.0);
         }
        }

   /// Muon Momentum at RE4 (Endcaps)

       id=24;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Momentum at RE4";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 30.0, 110.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }

   /// Loss of Muon Momentum in Iron (between RE1_layer1 and RE4)
       id=25;

    if (mh_RPC_tof.count(id) == 0) {
          ss<<"Loss of muon Momentum in Iron (endcap) ";
          mh_RPC_tof[id]  = new TH1F(ss.str().c_str(),"", 80, 0.0, 40.0);
          mh_RPC_tof[id]->SetXTitle(" Momentum (GeV/c)");
          mh_RPC_tof[id]->SetYTitle("Entries");
          ss.str("");       
    }
       if (station == 4 )
        {
         if ( touche4 == 0)
         {
          mome4=RPC[i]._pabs;
          touche4 = 1;
          mh_RPC_tof[24]->Fill(mome4,1.0);
          if (touche1 == 1 )
          {
            mh_RPC_tof[25]->Fill(mome1-mome4,1.0);
          }
         }
        }

  }  // End of Endcaps

   /// X-Local Coordinate vs Z-Local Coordinate

 //      id=15;
 //   if (mh_RPC_pos.count(id) == 0) {
 //         ss<<"Local x-coord. vs local z-coord of muon hit";
 //         mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -150., 150., 100, -0.1, 0.1 );
 //         mh_RPC_pos[id]->SetXTitle(" X-coord (cm)");
 //         mh_RPC_pos[id]->SetYTitle("Z-coord (cm)");
 //         ss.str("");       
 //   }
 //       mh_RPC_pos[id]->Fill(RPC[i]._locposx, RPC[i]._locposz,1.0);

   /// X-Local Coordinate vs Y-Local Coordinate

       id=16;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<"local x-coord. vs local y-coord of muon hit";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -150., 150., 100, -100., 100. );
          mh_RPC_pos[id]->SetXTitle(" X-local coord (cm)");
          mh_RPC_pos[id]->SetYTitle("Y-local coord (cm)");
          ss.str("");       
    }
        mh_RPC_pos[id]->Fill(RPC[i]._locposx, RPC[i]._locposy,1.0);

   /// X-Global Coordinate vs Z-Global Coordinate

    radius = RPC[i]._globposz* ( 1.+ exp(-2.*RPC[i]._globposeta) )
        / ( 1. - exp(-2.*RPC[i]._globposeta) ) ;

    costeta =( 1. - exp(-2.*RPC[i]._globposeta) )/( 1. + exp(-2.*RPC[i]._globposeta) );

    sinteta = 2. * exp(-RPC[i]._globposeta) /( 1. + exp(-2.*RPC[i]._globposeta) );
    xposglob = radius*sinteta*cos(RPC[i]._globposphi );
    yposglob = radius*sinteta*sin(RPC[i]._globposphi );

   /// Radius of hit

       id = 19;
       if (mh_RPC_enloss.count(id) == 0) {
          ss<<" radius of hit (barrel) ";
          mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 100, 0., 1200.);
          mh_RPC_enloss[id]->SetXTitle("Radius (cm)");
          mh_RPC_enloss[id]->SetYTitle("Entries");
          ss.str("");       
       }
     if (region == 0 )  mh_RPC_enloss[id]->Fill(radius,1.0);

         id = 30;
       if (mh_RPC_enloss.count(id) == 0) {
          ss<<" radius of hit (endcaps) ";
          mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 100, 0., 1300.);
          mh_RPC_enloss[id]->SetXTitle("Radius (cm)");
          mh_RPC_enloss[id]->SetYTitle("Entries");
          ss.str("");
       }
     if (region != 0 )  mh_RPC_enloss[id]->Fill(radius,1.0);

   /// Costheta of hit

       id = 20;
       if (mh_RPC_enloss.count(id) == 0) {
          ss<<" costheta of hit (barrel) ";
          mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 100, -1., 1.);
          mh_RPC_enloss[id]->SetXTitle(" cos(theta) ");
          mh_RPC_enloss[id]->SetYTitle("Entries");
          ss.str("");       
       }
     if (region == 0 )  mh_RPC_enloss[id]->Fill(costeta,1.0);

          id = 31;
       if (mh_RPC_enloss.count(id) == 0) {
          ss<<" costheta of hit (endcaps) ";
          mh_RPC_enloss[id]  = new TH1F(ss.str().c_str(),"", 100, -1., 1.);
          mh_RPC_enloss[id]->SetXTitle(" cos(theta) ");
          mh_RPC_enloss[id]->SetYTitle("Entries");
          ss.str("");
       }
     if (region != 0 )  mh_RPC_enloss[id]->Fill(costeta,1.0);

       id=21;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<" path followed by muons vs radius ";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 160, 0., 160., 100, 200., 1300. );
          mh_RPC_pos[id]->SetXTitle(" Path followed by muon");
          mh_RPC_pos[id]->SetYTitle("Radius (cm)");
          ss.str("");       
    }
        mh_RPC_pos[id]->Fill((float)path, radius,1.0);


       id=17;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global z-coord of muon hit (barrel) ";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -800., 800., 100, -800., 800. );
          mh_RPC_pos[id]->SetXTitle(" Global Z-coord (cm)");
          mh_RPC_pos[id]->SetYTitle("Global X-coord (cm)");
          ss.str("");       
    }
    if (region == 0 ) mh_RPC_pos[id]->Fill( RPC[i]._globposz,xposglob, 1.0);

      id=32;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global z-coord of muon hit (endcaps ) ";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -1100., 1100., 100, -800., 800. );
          mh_RPC_pos[id]->SetXTitle(" Global Z-coord (cm)");
          mh_RPC_pos[id]->SetYTitle("Global X-coord (cm)");
          ss.str("");
    }
    if (region != 0 )  mh_RPC_pos[id]->Fill( RPC[i]._globposz,xposglob, 1.0);


   /// X-Global Coordinate vs Y-Global Coordinate

       id=18;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global y-coord of muon hit (barrel) ";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -800., 800., 100, -800., 800. );
          mh_RPC_pos[id]->SetXTitle(" Global X-coord (cm)");
          mh_RPC_pos[id]->SetYTitle("Global Y-coord (cm)");
          ss.str("");       
    }
     if (region == 0 )   mh_RPC_pos[id]->Fill(xposglob, yposglob,1.0);

      id=33;
    if (mh_RPC_pos.count(id) == 0) {
          ss<<"Global x-coord. vs global y-coord of muon hit ( endcaps) ";
          mh_RPC_pos[id]  = new TH2F(ss.str().c_str(),"", 100, -800., 800., 100, -800., 800. );
          mh_RPC_pos[id]->SetXTitle(" Global X-coord (cm)");
          mh_RPC_pos[id]->SetYTitle("Global Y-coord (cm)");
          ss.str("");
    }
     if (region != 0 )   mh_RPC_pos[id]->Fill(xposglob, yposglob,1.0);

   }
  }
}
hist_file->Write();
delete hist_file;
}
