#include "TText.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void DoCompare_CSC( ){

 static const int NHisto = 38;

 TText* te = new TText();
 te->SetTextSize(0.1);
 
  gROOT->ProcessLine(".x HistoCompare.C");
  HistoCompare * myPV = new HistoCompare();

 char*  reffilename  = "${REFFILE}";//"./CSCSimHitsPlots_ref.root";
 char*  curfilename  = "${CURFILE}";//"./CSCSimHitsPlots.root";

 TFile * reffile = new TFile(reffilename);
 TFile * curfile = new TFile(curfilename);

 //1-Dimension Histogram
 char* label[NHisto];
 label[0] = "Number_of_all_CSC_hits";
 label[1] = "Number_of_muon_CSC_hits";

 TH1F* htemp1[NHisto];
 TH1F* htemp2[NHisto];

 char labelh[10];
 int nh =0;
 char histoName[40];
 char histoNametof[40];

 char title[50];
 TCanvas c1;


 for ( int k = 1; k <3 ; ++k) {
   for ( int i = 1; i<5 ; ++i)  {
     for ( int j = 1; j<5 ; ++j)  {

         if (i != 1 && j>2 ) continue;
         if (i == 4 && j>1 ) continue;
         nh++;
         cout << " nh " << nh << endl;
         
    // Energy Loss plots
         int idhisto = k*100+ i*10+j;         
         sprintf(labelh,"ME%d", idhisto) ;
         strcpy(histoName, labelh);
         strcat(histoName,"_energy_loss");
         label[2*nh] = histoName;
         cout << " label energy" << label[2*nh] << endl; 
 
         htemp1[2*nh]  = dynamic_cast<TH1F*>(reffile->Get(label[2*nh])); 
         htemp2[2*nh]  = dynamic_cast<TH1F*>(curfile->Get(label[2*nh]));
         if( htemp1[2*nh] == 0 ) std::cout << " reference histo is empty " << endl;
         if( htemp2[2*nh] == 0 ) std::cout << " current histo is empty " << endl;   
         if( htemp1[2*nh] == 0 || htemp2[2*nh] == 0) continue;

         htemp1[2*nh]->SetLineColor(2);
         htemp2[2*nh]->SetLineColor(4);
         htemp1[2*nh]->SetLineStyle(3);
         htemp2[2*nh]->SetLineStyle(5);
         TLegend leg(0.1, 0.15, 0.2, 0.25);
         leg.AddEntry(htemp1[2*nh], "Reference", "l");
         leg.AddEntry(htemp2[2*nh], "New ", "l");

         htemp1[2*nh]->Draw();

         htemp2[2*nh]->Draw("Same");
         leg.Draw();

         myPV->PVCompute(htemp1[2*nh],htemp2[2*nh], te);
         sprintf(title,"%s%s", label[2*nh],".eps");
         c1.Print(title);

    // ToF plots
         int idhisto_tof = idhisto+200;
         strcpy(histoNametof, labelh);
         strcat(histoNametof,"_tof");
         label[2*nh+1] = histoNametof;
         cout << " label tof" << label[2*nh+1] << endl;

         htemp1[2*nh+1]  = dynamic_cast<TH1F*>(reffile->Get(label[2*nh+1]));
         htemp2[2*nh+1]  = dynamic_cast<TH1F*>(curfile->Get(label[2*nh+1]));
         if( htemp1[2*nh+1] == 0 ) std::cout << " reference histo is empty " << endl;
         if( htemp2[2*nh+1] == 0 ) std::cout << " current histo is empty " << endl;
         if( htemp1[2*nh+1] == 0 || htemp2[2*nh+1] == 0) continue;

         htemp1[2*nh+1]->SetLineColor(2);
         htemp2[2*nh+1]->SetLineColor(4);
         htemp1[2*nh+1]->SetLineStyle(3);
         htemp2[2*nh+1]->SetLineStyle(5);
         TLegend leg(0.1, 0.15, 0.2, 0.25);
         leg.AddEntry(htemp1[2*nh+1], "Reference", "l");
         leg.AddEntry(htemp2[2*nh+1], "New ", "l");

         htemp1[2*nh+1]->Draw();

         htemp2[2*nh+1]->Draw("Same");
         leg.Draw();

         myPV->PVCompute(htemp1[2*nh+1],htemp2[2*nh+1], te);
         sprintf(title,"%s%s", label[2*nh+1],".eps");
         c1.Print(title);

     }
   }
 }


 for ( int i = 0; i< 2 ; i++ ) {

     htemp1[i]  = dynamic_cast<TH1F*>(reffile->Get(label[i]));
     htemp2[i]  = dynamic_cast<TH1F*>(curfile->Get(label[i]));
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
   
 }


}


