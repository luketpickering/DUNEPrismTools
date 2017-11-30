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

/*  gTree = (TTree*)fin->Get("gRooTracker");

  if(!gTree) std::cout << "Can't get gtree" << std::endl;
  #ifdef DEBUG
  else std::cout<<"Got gtree" << std::endl; 
  #endif*/

  nEntries = inTree->GetEntries();
  #ifdef DEBUG
  std::cout << nEntries << std::endl; 
  #endif

  #ifdef DEBUG
  std::cout << "Set branches" << std::endl;
  #endif

  eventTree = new TTree("events","Event Information");
  fout = new TFile(outFileName.c_str(),"RECREATE");
//  gOutTree = (TTree*)gTree->Clone("genie"); 

  SetBranches();

  dimension[0] = det[0]*100;//Convert to cm
  dimension[1] = det[1]*100;
  dimension[2] = det[2]*100;
  FV[0] = fid[0]*100;//Convert to cm
  FV[1] = fid[1]*100;
  FV[2] = fid[2]*100;
  shift = off*100;

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
    if(!(ie%100))std::cout<<"Entry " << ie<<std::endl;
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
/*    std::cout << "vtxX " << EvtVtx[0] << std::endl;
    std::cout << "vtxY " << fabs(EvtVtx[1]) << std::endl; 
    std::cout << "vtxZ " << fabs(EvtVtx[2]) << std::endl; */
    #endif 

    double wallX[2] = {
        (shift - dimension[0]*.5 + FV[0]),
        (shift + dimension[0]*.5 - FV[0])};
    double wallY = dimension[1]*.5 - FV[1];
    double wallZ = dimension[2]*.5 - FV[2];

    //First: check if in y & z FV
    if( fabs(EvtVtx[1]) > wallY ) continue;
    else if( fabs(EvtVtx[2]) > wallZ ) continue;

    //Check if in x 
    else if( EvtVtx[0] > wallX[1] || EvtVtx[0] < wallX[0]) continue;

    #ifdef DEBUG
//    std::cout << "Inside FV" << std::endl;
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
    flagExitBack = false;
    flagExitFront = false;
    flagExitY = false;
    flagExitXHigh = false;
    flagExitXLow = false;
    flagMuContained = false;
    flagNoEHadOut = false;

    muExitingPX = 0.;
    muExitingPY = 0.;
    muExitingPZ = 0.;

    eHadOutDep = 0.;
    eHadInDep = 0.;
    eHadTotalDep = 0.;

    eMuDep = 0.;
    eMuTotalDep = 0.;

    std::vector<int>::iterator itPi0;
    std::map<int,int>::iterator itGamma;

    trackPi0.clear();
    trackGamma.clear();
    trackE.clear();
    trackEIn.clear();
    trackEOut.clear();

    ///////Start of step loop
    for(int i = 0; i < nstep; ++i){

      //Look for muon track
      if(PID[i] == 13 && !muEnd){
        if(!found) {
//          std::cout<<"found start of muon track"<<std::endl;
          found = true;          
        }                
        
        eMuDep += edep[i];
        eMuTotalDep += edep[i];

        //check if it hits the walls
        if(xe[i] < wallX[0] - FV[0]){
          flagExitXLow = true; 
          muExitingPX = pxe[i];
          muExitingPY = pye[i];
          muExitingPZ = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted XLow "<< 
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
        else if(xe[i] > wallX[1] + FV[0]){ 
          flagExitXHigh = true; 
          muExitingPX = pxe[i];
          muExitingPY = pye[i];
          muExitingPZ = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted right "<< 
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }

        if(abs(ye[i]) > wallY + FV[1]){ 
          flagExitY = true; 
          muExitingPX = pxe[i];
          muExitingPY = pye[i];
          muExitingPZ = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted Y "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }

        if(ze[i] > wallZ + FV[2]){ 
          flagExitBack = true; 
          muExitingPX = pxe[i];
          muExitingPY = pye[i];
          muExitingPZ = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted back "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
        else if(ze[i] < -wallZ - FV[2]){ 
          flagExitFront = true; 
          muExitingPX = pxe[i];
          muExitingPY = pye[i];
          muExitingPZ = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted front "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
      }

      else if(PID[i] == 13 && muEnd){//Energy dep outside detector
        eMuTotalDep += edep[i];
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
            
            eHadOutDep += edep[i]; 
          }

          else eHadInDep += edep[i];
        }
      }

      //Look for Pi0. Add in track id if not in vector
      else if(PID[i] == 111){
     //   std::cout << "Found pi0" << std::endl;
        itPi0 = std::find( trackPi0.begin(), trackPi0.end(), track[i] );
        if(itPi0 == trackPi0.end() ){
       //   std::cout << "\tNew pi0 "<< track[i]<<std::endl;
          trackPi0.push_back(track[i]);
        }
      }

      //Look for gamma. Add in track and parid if not in vector
      else if(PID[i] == 22){
    //    std::cout<<"Found gamma "<<parid[i] << std::endl;
        if(trackGamma.find(track[i]) == trackGamma.end() ){
   //       std::cout << "\t new gamma " << track[i] << " "  << parid[i]<<std::endl;
          trackGamma[track[i]] = parid[i];
        }
      }

      //Add in all edeps and parids from electrons
      else if(abs(PID[i]) == 11){        
        trackE.push_back( std::make_pair(parid[i],edep[i]) );

        //Inside walls. Excuse terrible naming
        if(xe[i] > (wallX[0] - FV[0]) && xe[i] < (wallX[1] + FV[0])&&
          abs(ye[i]) < (wallY + FV[0]) && abs(ze[i]) < (wallZ + FV[0])
          ){
        
          //Between wall and FV
          if(xe[i] < wallX[0] || xe[i] > wallX[1] ||
            abs(ye[i]) > wallY || abs(ze[i]) > wallZ){
            
            trackEOut.push_back( std::make_pair(parid[i],edep[i]) ); 
          }

          else trackEIn.push_back( std::make_pair(parid[i],edep[i]) );
        }
      }    
    }

    flagMuContained = !(flagExitFront || flagExitBack ||
                         flagExitXLow || flagExitXHigh ||
                         flagExitY);
        
    if(trackE.size() == 0)std::cout<<"found no e+-"<<std::endl;
    
    ePi0TotalDep = 0.; 
    ePi0InDep = 0.;
    ePi0OutDep = 0.;

    //Iterate over e+- edeps, check if parid in gamma map
    //if so: check if gamma parid in pi0 vector
    for(int i = 0; i < trackE.size(); ++i){
      int eParent = trackE.at(i).first;
      if(trackGamma.find(eParent) != trackGamma.end()){ // if in
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
 //         std::cout << "FOUND PI0->G->E" <<std::endl;
          ePi0TotalDep += trackE.at(i).second; //add in edep
        }
      }
    }

    for(int i = 0; i < trackEIn.size(); ++i){
      int eParent = trackEIn.at(i).first;
      if(trackGamma.find(eParent) != trackGamma.end()){ // if in
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
 //         std::cout << "FOUND PI0->G->E" <<std::endl;
          ePi0InDep += trackEIn.at(i).second; //add in edep
        }
      }
    }

    for(int i = 0; i < trackEOut.size(); ++i){
      int eParent = trackEOut.at(i).first;
      if(trackGamma.find(eParent) != trackGamma.end()){ // if in
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
 //         std::cout << "FOUND PI0->G->E" <<std::endl;
          ePi0OutDep += trackEOut.at(i).second; //add in edep
        }
      }
    }

    flagNoEHadOut = (ePi0OutDep == 0 && ( ePi0TotalDep == ePi0InDep )
                     && eHadOutDep == 0 
                     && ( eHadTotalDep == eHadInDep ) );
    ////////////////////////////////////

    eHadTrueCharged = 0.;
    eHadTrueTotal = 0.;
    eMuTrue = 0.;

    nMu = 0;
    nPi0 = 0;
    nPiC = 0;
    nProton = 0;
    nNeutron = 0;

    //Start of particle loop
//    std::cout << "Starting true FS particle loop "<< ni << std::endl;

    for(int ip = 0; ip < ni; ++ip){
      //std::cout << ip << std::endl;
      if(PIDi[ip] == 13){
        eMuTrue += ekini[ip] + mi[ip];
        nMu++;
      }
      else if(PIDi[ip] == 111){
        eHadTrueCharged += ekini[ip] + mi[ip];
        nPi0++;
      }
      else if(abs(PIDi[ip]) == 211 ){
        eHadTrueCharged += ekini[ip] + mi[ip];
        nPiC++;
      }
      else if(PIDi[ip] == 2212){
        eHadTrueCharged += ekini[ip];
        nProton++;
      }
      else if(PIDi[ip] == 2112){
        eHadTrueTotal += ekini[ip];
        nNeutron++;
      }
    }
    eHadTrueTotal += eHadTrueCharged; 
    ///////////////////////////////////

    eventNum = ev; 
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

  inTree->SetBranchAddress("ev",&ev);
  inTree->SetBranchAddress("ekina",&Enu);
  inTree->SetBranchAddress("xa",&EvtVtx[0]);
  inTree->SetBranchAddress("ya",&EvtVtx[1]);
  inTree->SetBranchAddress("za",&EvtVtx[2]);
  inTree->SetBranchAddress("pida",&nuPID);
  
  inTree->SetBranchAddress("pid",&PID);
  inTree->SetBranchAddress("tid",&track);
  inTree->SetBranchAddress("parid",&parid);
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
 
  eventTree->Branch("eventNum",&eventNum);
  eventTree->Branch("Enu",&Enu,"Enu/D");
  eventTree->Branch("vtx_X",&EvtVtx[0],"vtx_X/D");
  eventTree->Branch("vtx_Y",&EvtVtx[1],"vtx_Y/D");
  eventTree->Branch("vtx_Z",&EvtVtx[2],"vtx_Z/D");

  eventTree->Branch("flagExitBack",&flagExitBack);
  eventTree->Branch("flagExitFront",&flagExitFront);
  eventTree->Branch("flagExitY",&flagExitY);
  eventTree->Branch("flagExitXLow",&flagExitXLow);
  eventTree->Branch("flagExitXHigh",&flagExitXHigh);
  eventTree->Branch("flagNoEHadOut",&flagNoEHadOut);
  eventTree->Branch("flagMuContained",&flagMuContained);

  eventTree->Branch("muExitingPX",&muExitingPX,"muExitingPX/D");
  eventTree->Branch("muExitingPY",&muExitingPY,"muExitingPY/D");
  eventTree->Branch("muExitingPZ",&muExitingPZ,"muExitingPZ/D");
  eventTree->Branch("eHadOutDep",&eHadOutDep,"eHadOutDep/D");
  eventTree->Branch("eHadInDep",&eHadInDep,"eHadInDep/D");
  eventTree->Branch("eHadTotalDep",&eHadTotalDep,"eHadTotalDep/D");
  eventTree->Branch("eMuDep",&eMuDep,"eMuDep/D");
  eventTree->Branch("eMuTotalDep",&eMuTotalDep,"eMuTotalDep/D");
  eventTree->Branch("ePi0TotalDep",&ePi0TotalDep,"ePi0TotalDep/D");
  eventTree->Branch("ePi0InDep",&ePi0InDep,"ePi0InDep/D");
  eventTree->Branch("ePi0OutDep",&ePi0OutDep,"ePi0OutDep/D");
  eventTree->Branch("eHadTrueCharged",&eHadTrueCharged,"eHadTrueCharged/D");
  eventTree->Branch("eHadTrueTotal",&eHadTrueTotal,"eHadTrueTotal/D");
  eventTree->Branch("eMuTrue",&eMuTrue,"eMuTrue/D");
  eventTree->Branch("nMu",&nMu,"nMu/I");
  eventTree->Branch("nPi0",&nPi0,"nPi0/I");
  eventTree->Branch("nPiC",&nPiC,"nPiC/I");
  eventTree->Branch("nProton",&nProton,"nProton/I");
  eventTree->Branch("nNeutron",&nNeutron,"nNeutron/I");

}

void DunePrismAnalyzer::Finalize(){
  fout->cd();
  eventTree->Write();
  //gOutTree->Write();
  
  fout->Close();
  fin->Close();

}
