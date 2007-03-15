#include "TText.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

void DoCompare_RPC( ){

 static const int NHisto = 33;

 TText* te = new TText();
 te->SetTextSize(0.1);
 
  gROOT->ProcessLine(".x HistoCompare.C");
  HistoCompare * myPV = new HistoCompare();

 char*  reffilename  = "${REFFILE}";//"./RPCSimHitsPlots_ref.root";
 char*  curfilename  = "${CURFILE}";//"./RPCSimHitsPlots.root";

 TFile * reffile = new TFile(reffilename);
 TFile * curfile = new TFile(curfilename);

 //1-Dimension Histogram
 char* label[NHisto];
 label[0] = "Number_of_all_RPC_hits";
 label[1] = "Number_of_muon_RPC_hits";
 label[2] = "Region_occupancy";
 label[3] = "Ring_occupancy_barrel";
 label[4] = "Ring_occupancy_endcaps";
 label[5] = "Station_occupancy_barrel";
 label[6] = "Station_occupancy_endcaps";
 label[7] = "Sector_occupancy_barrel";
 label[8] = "Sector_occupancy_endcaps";
 label[9] = "Layer_occupancy_barrel";
 label[10] = "Layer_occupancy_endcaps";
 label[11] = "Subsector_occupancy_barrel";
 label[12] = "Subsector_occupancy_endcaps";
 label[13] = "Roll_occupancy_barrel";
 label[14] = "Roll_occupancy_endcaps";
 label[15] = "RPC_energy_loss_barrel";
 label[16] = "RPC_energy_loss_endcaps";
 label[17] = "path_followed_by_muon";
 label[18] = "Momentum_at_RB1";
 label[19] = "Momentum_at_RB4";
 label[20] = "Loss_of_muon_Momentum_in_Iron_barrel" ;
 label[21] = "Momentum_at_RE1";                                
 label[22] = "Momentum_at_RE4";                                 
 label[23] = "Loss_of_muon_Momentum_in_Iron_endcap";
 label[24] = "radius_of_hit_barrel";
 label[25] = "radius_of_hit_endcaps";
 label[26] = "costheta_of_hit_barrel";
 label[27] = "costheta_of_hit_endcaps";
 label[28] = "local_x-coord_vs_local_y-coord_of_muon_hit";
 label[29] = "Global_z-coord_vs_global_x-coord_of_muon_hit_barrel";
 label[30] = "Global_x-coord_vs_global_y-coord_of_muon_hit_barrel";
 label[31] = "Global_z-coord_vs_global_x-coord_of_muon_hit_endcaps";
 label[32] = "Global_x-coord_vs_global_y-coord_of_muon_hit_endcaps";



 TH1F* htemp1[NHisto];
 TH1F* htemp2[NHisto];

 for ( int i = 0; i< NHisto ; i++ ) {
   char title[70];
   TCanvas c1;

   if ( i<28 ) 
   {
     htemp1[i]  = dynamic_cast<TH1F*>(reffile->Get(label[i]));
     htemp2[i]  = dynamic_cast<TH1F*>(curfile->Get(label[i]));
     if( htemp1[i] == 0 ) std::cout << " reference histo is empty " << endl;
     if( htemp2[i] == 0 ) std::cout << " current histo is empty " << endl;
     if( htemp1[i] == 0 || htemp2[i] == 0) continue;

     htemp1[i]->SetLineColor(2);
     htemp2[i]->SetLineColor(4);
     htemp1[i]->SetLineStyle(3);
     htemp2[i]->SetLineStyle(5);
     TLegend leg(0.1, 0.15, 0.2, 0.25);
     leg.AddEntry(htemp1[i], "Reference", "l");
     leg.AddEntry(htemp2[i], "New ", "l");

     htemp1[i]->Draw();

     htemp2[i]->Draw("Same"); 
     leg.Draw();
     myPV->PVCompute(htemp1[i],htemp2[i], te);
     sprintf(title,"%s%s", label[i],".eps");
     c1.Print(title);
   
   } else {

     htemp1[i]  = dynamic_cast<TH2F*>(reffile->Get(label[i]));
     htemp2[i]  = dynamic_cast<TH2F*>(curfile->Get(label[i]));

     htemp1[i]->SetMarkerStyle(21);
     htemp2[i]->SetMarkerStyle(22);
     htemp1[i]->SetMarkerColor(2);
     htemp2[i]->SetMarkerColor(4);
     htemp1[i]->SetMarkerSize(0.3);
     htemp2[i]->SetMarkerSize(0.3);
 
     c1.Divide(1,2);
     c1.cd(1);
     htemp1[i]->Draw();
     leg.Draw();
     
     c1.cd(2);
     htemp2[i]->Draw();
     leg.Draw();

     sprintf(title,"%s%s", label[i],".eps");  
     c1.Print(title);
  } 
 }


}


