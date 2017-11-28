#include "DunePrismAnalyzer.h"

#ifndef DEBUG
 #define DEBUG
#endif
int main(int argc, char* argv[]){
  parse_args(argc, argv);

  DunePrismAnalyzer dpa = DunePrismAnalyzer(inFileName, outFileName, detector, fiducialGap, offset);
  dpa.Analyze();
  dpa.Finalize();
  return 0;
}

DunePrismAnalyzer::DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off){
  fin = new TFile(inFileName.c_str(),"READ");
  inTree = (TTree*)fin->Get("argon"); 
  gTree = (TTree*)fin->Get("gRooTracker");
  nEntries = inTree->GetEntries();
 
  Enu = 0.;
  EvtVtx[0] = 0.;
  EvtVtx[1] = 0.;
  EvtVtx[2] = 0.;
  nuPID = 0;
  inTree->SetBranchAddress("ekina",&Enu);
  inTree->SetBranchAddress("xa",&EvtVtx[0]);
  inTree->SetBranchAddress("ya",&EvtVtx[1]);
  inTree->SetBranchAddress("za",&EvtVtx[2]);
  inTree->SetBranchAddress("pida",&nuPID);

  fout = new TFile(outFileName.c_str(),"RECREATE");
  eventTree = new TTree("events","Event Information");

  eventTree->Branch("Enu",&Enu,"Enu/D");
  eventTree->Branch("vtx_X",&EvtVtx[0],"vtx_X/D");
  eventTree->Branch("vtx_Y",&EvtVtx[1],"vtx_Y/D");
  eventTree->Branch("vtx_Z",&EvtVtx[2],"vtx_Z/D");

  gOutTree = (TTree*)gTree->Clone("genie"); 

  detector[0] = det[0]*1.E3;//Convert to mm
  detector[1] = det[1]*1.E3;
  detector[2] = det[2]*1.E3;
  fiducialGap[0] = fid[0]*1.E3;//Convert to mm
  fiducialGap[1] = fid[1]*1.E3;
  fiducialGap[2] = fid[2]*1.E3;
  offset = off*1.E3;

}

void DunePrismAnalyzer::Analyze(){
  for (int i = 0; i < nEntries; ++i){
    inTree->GetEntry(i);
    if(nuPID != 14) continue;
    #ifdef DEBUG
    std::cout << i << std::endl;
    std::cout<<Enu<<std::endl;
    std::cout<<nuPID<<std::endl;
    #endif    
 
    //To Do: make a vector of trees
    //that correspond to different 
    //stops off axis
    
    //First: check if in y & z FV
    if( abs(EvtVtx[1]) > (detector[1]*.5 - fiducialGap[1]) ) continue;
    if( abs(EvtVtx[2]) > (detector[2]*.5 - fiducialGap[2]) ) continue;

    if( EvtVtx[0] > (offset + detector[0]/2 - fiducialGap[0]) || EvtVtx[0] < (offset - detector[0]/2 + fiducialGap[0])) continue;

    //Check if in 
    eventTree->Fill();
  }
}

void DunePrismAnalyzer::Finalize(){
  fout->cd();
  eventTree->Write();
  gOutTree->Write();
  
  fout->Close();
  fin->Close();

}
