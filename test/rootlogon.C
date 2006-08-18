{
//gSystem->Load("libCintex.so");   // for CMSSW_0_8_0_pre4
//Cintex::Enable();                 // 

gSystem->Load("libFWCoreFWLite");
AutoLibraryLoader::enable();

#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>

TFile* hist_file; 

gROOT->ProcessLine(".L mymap.C+");
}
