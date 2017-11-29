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
  #ifdef DEBUG
  std::cout << inFileName << std::endl;
  std::cout << outFileName << std::endl;
  std::cout << det[0] << " " << det[1] << " " << det[2] << std::endl;
  std::cout << fid[0] << " " << fid[1] << " " << fid[2] << std::endl;
  std::cout << off << std::endl;
  #endif

  fin = new TFile(inFileName.c_str(),"READ");
  inTree = (TTree*)fin->Get("argon"); 

  if(!inTree) std::cout << "Can't get argon tree" << std::endl;
  #ifdef DEBUG
  else std::cout<<"Got argon tree" << std::endl;
  #endif

  gTree = (TTree*)fin->Get("gRooTracker");

  if(!gTree) std::cout << "Can't get gtree" << std::endl;
  #ifdef DEBUG
  else std::cout<<"Got gtree" << std::endl; 
  #endif

  nEntries = inTree->GetEntries();
  #ifdef DEBUG
  std::cout << nEntries << std::endl; 
  #endif

  #ifdef DEBUG
  std::cout << "Set branches" << std::endl;
  #endif

  eventTree = new TTree("events","Event Information");
  fout = new TFile(outFileName.c_str(),"RECREATE");
  gOutTree = (TTree*)gTree->Clone("genie"); 

  SetBranches();

  dimension[0] = det[0]*1.E3;//Convert to mm
  dimension[1] = det[1]*1.E3;
  dimension[2] = det[2]*1.E3;
  FV[0] = fid[0]*1.E3;//Convert to mm
  FV[1] = fid[1]*1.E3;
  FV[2] = fid[2]*1.E3;
  shift = off*1.E3;

  #ifdef DEBUG
  std::cout << dimension[0] << std::endl;
  std::cout << dimension[1] << std::endl;
  std::cout << dimension[2] << std::endl;
  std::cout << FV[0] << std::endl; 
  std::cout << FV[1] << std::endl;
  std::cout << FV[2] << std::endl;
  std::cout << shift << std::endl;
  #endif
}
 
void DunePrismAnalyzer::Analyze(){
  for (int ie = 0; ie < nEntries; ++ie){
    inTree->GetEntry(ie);
    if(nuPID != 14) continue;
/*    #ifdef DEBUG
    std::cout << i << std::endl;
    std::cout<<Enu<<std::endl;
    std::cout<<nuPID<<std::endl;
    #endif */
 
    //To Do: make a vector of trees
    //that correspond to different 
    //stops off axis
    
    //////////////Vertex in Detector?
    #ifdef DEBUG
    std::cout << "vtxX " << EvtVtx[0] << std::endl;
    std::cout << "vtxY " << fabs(EvtVtx[1]) << std::endl; 
    std::cout << "vtxZ " << fabs(EvtVtx[2]) << std::endl; 
    #endif 

    double wallX[2] = {
        (shift - dimension[0]/2 + FV[0]),
        (shift + dimension[0]/2 - FV[0])};
    double wallY = dimension[1]*.5 - FV[1];
    double wallZ = dimension[2]*.5 - FV[2];

    //First: check if in y & z FV
    if( fabs(EvtVtx[1]) > wallY ) continue;
    else if( fabs(EvtVtx[2]) > wallZ ) continue;

    //Check if in x 
    else if( EvtVtx[0] > wallX[1] || EvtVtx[0] < wallX[0]) continue;

    #ifdef DEBUG
    std::cout << "Inside FV" << std::endl;
    #endif
    /////////////////////////////////////
    
    mu_x.clear();
    mu_y.clear();
    mu_z.clear();
    mu_px.clear();
    mu_py.clear();
    mu_pz.clear();
    
    bool found = false;
    bool muEnd = false;
    flagExitBack = 0;
    flagExitFront = 0;
    flagExitY = 0;
    flagExitRight = 0;
    flagExitLeft = 0;

    muExitingP[0] = 0.;
    muExitingP[1] = 0.;
    muExitingP[2] = 0.;

    eHadOut = 0.;
    eHadIn = 0.;
    eHadTotalDep = 0.;
    ///////Start of step loop
    for(int i = 0; i < nstep; ++i){

      //Look for muon track
      if(PID[i] == 13 && !muEnd){
        if(!found) {
          std::cout<<"found start of muon track"<<std::endl;
          found = true;          
        }                

        //check if it hits the walls
        if(xe[i] < wallX[0] - FV[0]){
          flagExitLeft = 1; 
          muExitingP[0] = pxe[i];
          muExitingP[1] = pye[i];
          muExitingP[2] = pze[i];
          std::cout<<"Muon exitted left"<< std::endl;
          muEnd = true;
        }
        else if(xe[i] > wallX[1] + FV[0]){ 
          flagExitRight = 1; 
          muExitingP[0] = pxe[i];
          muExitingP[1] = pye[i];
          muExitingP[2] = pze[i];
          std::cout<<"Muon exitted right"<< std::endl;
          muEnd = true;
        }
        else if(ye[i] > wallY + FV[1]){ 
          flagExitY = 1; 
          muExitingP[0] = pxe[i];
          muExitingP[1] = pye[i];
          muExitingP[2] = pze[i];
          std::cout<<"Muon exitted Y"<< std::endl;
          muEnd = true;
        }
        else if(ze[i] > wallZ + FV[2]){ 
          flagExitBack = 1; 
          muExitingP[0] = pxe[i];
          muExitingP[1] = pye[i];
          muExitingP[2] = pze[i];
          std::cout<<"Muon exitted back"<< std::endl;
          muEnd = true;
        }
        else if(ze[i] < -wallZ - FV[2]){ 
          flagExitFront = 1; 
          muExitingP[0] = pxe[i];
          muExitingP[1] = pye[i];
          muExitingP[2] = pze[i];
          std::cout<<"Muon exitted front"<< std::endl;
          muEnd = true;
        }
      }

      //Look for hadronic energy deposits
      //right now: protons and charged pions
      //need to add in em showers from pi0
      else if(PID[i] == 2212 || abs(PID[i]) == 211){
        //std::cout << "checking hadron edeps" << std::endl;        
        eHadTotalDep += edep[i];

        //Inside walls. Excuse terrible naming
        if(xe[i] > (wallX[0] - FV[0]) && xe[i] < (wallX[1] + FV[0])&&
          abs(ye[i]) < (wallY + FV[0]) && abs(ze[i]) < (wallZ + FV[0])
          ){
        
          //Between wall and FV
          if(xe[i] < wallX[0] || xe[i] > wallX[1] ||
            abs(ye[i]) > wallY || abs(ze[i]) > wallZ){
            
            eHadOut += edep[i]; 
          }

          else eHadIn += edep[i];
        }
      }
    }
    ////////////////////////////////////

    eHadTrueCharged = 0.;
    eHadTrueTotal = 0.;
    eMu = 0.;

    //Start of particle loop
    std::cout << "Starting true FS particle loop "<< ni << std::endl;

    for(int ip = 0; ip < ni; ++ip){
      //std::cout << ip << std::endl;
      if(PIDi[ip] == 13) eMu += ekini[ip] + mi[ip];
      else if(PIDi[ip] == 111 || abs(PIDi[ip]) == 211 ) eHadTrueCharged += ekini[ip] + mi[ip]; 
      else if(PIDi[ip] == 2212) eHadTrueCharged += ekini[ip];
      else if(PIDi[ip] == 2112) eHadTrueTotal += ekini[ip];
    }
    eHadTrueTotal += eHadTrueCharged; 
    ///////////////////////////////////
    eventTree->Fill();
  }
}

void DunePrismAnalyzer::SetBranches(){

  inTree->SetBranchAddress("ni",&ni);
  inTree->SetBranchAddress("pidi",&PIDi);
  inTree->SetBranchAddress("pxi",&pxi);
  inTree->SetBranchAddress("pyi",&pyi);
  inTree->SetBranchAddress("pzi",&pzi);
  inTree->SetBranchAddress("ekini",&ekini);
  inTree->SetBranchAddress("mi",&mi);

  inTree->SetBranchAddress("ekina",&Enu);
  inTree->SetBranchAddress("xa",&EvtVtx[0]);
  inTree->SetBranchAddress("ya",&EvtVtx[1]);
  inTree->SetBranchAddress("za",&EvtVtx[2]);
  inTree->SetBranchAddress("pida",&nuPID);
  
  inTree->SetBranchAddress("pid",&PID);
  inTree->SetBranchAddress("nstep",&nstep);
  inTree->SetBranchAddress("xs",&xs);
  inTree->SetBranchAddress("xe",&xe);
  inTree->SetBranchAddress("ys",&ys);
  inTree->SetBranchAddress("ye",&ye);
  inTree->SetBranchAddress("zs",&zs);
  inTree->SetBranchAddress("ze",&ze);
  inTree->SetBranchAddress("pxs",&pxs);
  inTree->SetBranchAddress("pxe",&pxe);
  inTree->SetBranchAddress("pys",&pys);
  inTree->SetBranchAddress("pye",&pye);
  inTree->SetBranchAddress("pzs",&pzs);
  inTree->SetBranchAddress("pze",&pze);
  inTree->SetBranchAddress("ekin",&ekin);
  inTree->SetBranchAddress("edep",&edep);
 
  eventTree->Branch("Enu",&Enu,"Enu/D");
  eventTree->Branch("vtx_X",&EvtVtx[0],"vtx_X/D");
  eventTree->Branch("vtx_Y",&EvtVtx[1],"vtx_Y/D");
  eventTree->Branch("vtx_Z",&EvtVtx[2],"vtx_Z/D");

  eventTree->Branch("flagExitBack",&flagExitBack);
  eventTree->Branch("flagExitFront",&flagExitFront);
  eventTree->Branch("flagExitY",&flagExitY);
  eventTree->Branch("flagExitLeft",&flagExitLeft);
  eventTree->Branch("flagExitRight",&flagExitRight);

  eventTree->Branch("muExitingP",&muExitingP,"muExitingP[3]/D");
  eventTree->Branch("eHadOut",&eHadOut,"eHadOut/D");
  eventTree->Branch("eHadIn",&eHadIn,"eHadIn/D");
  eventTree->Branch("eHadTotalDep",&eHadTotalDep,"eHadTotalDep/D");
  eventTree->Branch("eHadTrueCharged",&eHadTrueCharged,"eHadTrueCharged/D");
  eventTree->Branch("eHadTrueTotal",&eHadTrueTotal,"eHadTrueTotal/D");
  eventTree->Branch("eMu",&eMu,"eMu/D");

}

void DunePrismAnalyzer::Finalize(){
  fout->cd();
  eventTree->Write();
  //gOutTree->Write();
  
  fout->Close();
  fin->Close();

}
