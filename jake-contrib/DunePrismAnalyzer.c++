#include "DunePrismAnalyzer.h"

#ifndef DEBUG
 #define DEBUG
#endif



int main(int argc, char* argv[]){
  parse_args_xml(argc, argv);



//  DunePrismAnalyzer dpa = DunePrismAnalyzer(inFileName, outFileName, detector, fiducialGap, offset);
  DunePrismAnalyzer dpa = DunePrismAnalyzer(inFileName, outFileName, detStops);
  std::cout << "Checking dpa size: " << dpa.stopsTrees.size()<< std::endl;
  dpa.AnalyzeStops();
  dpa.FinalizeStops();
  return 0;
}

/*DunePrismAnalyzer::DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off){
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

  nEntries = inTree->GetEntries();
  #ifdef DEBUG
  std::cout << nEntries << std::endl; 
  #endif

  #ifdef DEBUG
  std::cout << "Set branches" << std::endl;
  #endif

  eventTree = new TTree("events","Event Information");
  fout = new TFile(outFileName.c_str(),"RECREATE");

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
}*/

DunePrismAnalyzer::DunePrismAnalyzer(std::string inFileName, std::string outFileName, std::vector<DetectorStop> detStops){
 #ifdef DEBUG
  std::cout << inFileName << std::endl;
  std::cout << outFileName << std::endl;
  std::cout << "Nstops: " << detStops.size() << std::endl;
  #endif

  fin = new TFile(inFileName.c_str(),"READ");
  inTree = (TTree*)fin->Get("argon"); 

  if(!inTree) std::cout << "Can't get argon tree" << std::endl;
  #ifdef DEBUG
  else std::cout<<"Got argon tree" << std::endl;
  #endif

  nEntries = inTree->GetEntries();
  #ifdef DEBUG
  std::cout << nEntries << std::endl; 
  #endif

  #ifdef DEBUG
  std::cout << "Set branches" << std::endl;
  #endif
  
  fout = new TFile(outFileName.c_str(),"RECREATE");

  nStops = detStops.size();
  std::cout << detStops.at(0).detectorSizeX << std::endl;

  for(int i = 0; i < nStops; ++i){
    InitVarsStop();
  }

  for(int i = 0; i < nStops; ++i){


    std::string treeName = "events_";
/*    std::string temp;

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << detStops.at(i).detectorSizeX;
    temp = stream.str();
    treeName += temp;

    treeName += "x";
    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) <<  detStops.at(i).detectorSizeY;
    temp = stream.str();
    treeName += temp;

    treeName += "x";

    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) << detStops.at(i).detectorSizeZ;
    temp = stream.str();
    treeName += temp;

    treeName += "_";

    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) << detStops.at(i).fiducialGapX;
    temp = stream.str();
    treeName += temp;

    treeName += "x";

    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) << detStops.at(i).fiducialGapY;
    temp = stream.str();
    treeName += temp;

    treeName += "x";

    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) << detStops.at(i).fiducialGapZ;
    temp = stream.str();
    treeName += temp;

    treeName += "_";

    stream.str(std::string());
    stream << std::fixed << std::setprecision(2) << detStops.at(i).shift;
    temp = stream.str();
    treeName += temp;
*/

    treeName += std::to_string(int(detStops.at(i).detectorSizeX*100));    
    treeName += "x";
    treeName += std::to_string(int(detStops.at(i).detectorSizeY*100));    
    treeName += "x";
    treeName += std::to_string(int(detStops.at(i).detectorSizeZ*100));    
    treeName += "_";

    treeName += std::to_string(int(detStops.at(i).fiducialGapX*100));    
    treeName += "x";
    treeName += std::to_string(int(detStops.at(i).fiducialGapY*100));    
    treeName += "x";
    treeName += std::to_string(int(detStops.at(i).fiducialGapZ*100));    
    treeName += "_";

    treeName += std::to_string(int(detStops.at(i).shift*100));

    std::cout << treeName << std::endl;

    TTree * sTree = new TTree(treeName.c_str(),"Event Information");
    stopsTrees.push_back(sTree);
    SetOutBranches(i);

    std::array<double,3> stop_dim = {detStops.at(i).detectorSizeX*100.,
                          detStops.at(i).detectorSizeY*100.,
                          detStops.at(i).detectorSizeZ*100.};

    std::array<double, 3> stop_FV = {detStops.at(i).fiducialGapX*100.,
                          detStops.at(i).fiducialGapY*100.,
                          detStops.at(i).fiducialGapZ*100.};

    double stop_shift = detStops.at(i).shift*100.;

    dimension.push_back(stop_dim); 
    FV.push_back(stop_FV);
    shift.push_back(stop_shift);

  }

  SetInBranches();
/*
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
*/
}

/*void DunePrismAnalyzer::Analyze(){
  for (int ie = 0; ie < nEntries; ++ie){
    if(!(ie%100))std::cout<<"Entry " << ie<<std::endl;
    inTree->GetEntry(ie);
    if(nuPID != 14) continue;
//    #ifdef DEBUG
//    std::cout << i << std::endl;
//    std::cout<<Enu<<std::endl;
//    std::cout<<nuPID<<std::endl;
//    #endif 
 
    //To Do: make a vector of trees
    //that correspond to different 
    //stops off axis
    
    //////////////Vertex in Detector?
    #ifdef DEBUG
//    std::cout << "vtxX " << EvtVtx[0] << std::endl;
//    std::cout << "vtxY " << fabs(EvtVtx[1]) << std::endl; 
//    std::cout << "vtxZ " << fabs(EvtVtx[2]) << std::endl; 
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
    std::vector<int>::iterator itPiC;
    std::vector<int>::iterator itProton;
    std::vector<int>::iterator itMu;

    ///New
    trackMu.clear();
    trackPiC.clear();
    trackProton.clear();
    
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
          itMu = std::find( trackMu.begin(), trackMu.end(), track[i] );
          if(itMu == trackMu.end() ){
            trackMu.push_back(track[i]);
          }
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
        if(PID[i] == 2212){
          itProton = std::find( trackProton.begin(), trackProton.end(), track[i] );
          if(itProton == trackProton.end() ){
            trackProton.push_back(track[i]);
          }
        }
        else{
          itPiC = std::find( trackPiC.begin(), trackPiC.end(), track[i] );
          if(itPiC == trackPiC.end() ){
            trackPiC.push_back(track[i]);
          }
        }
      }

      //Look for Pi0. Add in track id if not in vector
      else if(PID[i] == 111){
        itPi0 = std::find( trackPi0.begin(), trackPi0.end(), track[i] );
        if(itPi0 == trackPi0.end() ){
          trackPi0.push_back(track[i]);
        }
      }

      //Look for gamma. Add in track and parid if not in map
      else if(PID[i] == 22){
        if(trackGamma.find(track[i]) == trackGamma.end() ){
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

    ePiCEMTotalDep = 0.; 
    ePiCEMInDep = 0.;
    ePiCEMOutDep = 0.;

    eProtonEMTotalDep = 0.; 
    eProtonEMInDep = 0.;
    eProtonEMOutDep = 0.;

    eMuEMTotalDep = 0.; 
    eMuEMInDep = 0.;
    eMuEMOutDep = 0.;

    //Iterate over e+- edeps, check if parid in gamma map
    //if so: check if gamma parid in pi0 vector
    for(int i = 0; i < trackE.size(); ++i){
      int eParent = trackE.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if parent = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0TotalDep += trackE.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMTotalDep += trackE.at(i).second;
//          std::cout << "Found pic->gamma->e"<<std::endl;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMTotalDep += trackE.at(i).second;
//          std::cout << "Found proton->gamma->e"<<std::endl;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMTotalDep += trackE.at(i).second;
 //         std::cout << "Found mu->gamma->e"<<std::endl;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMTotalDep += trackE.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMTotalDep += trackE.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMTotalDep += trackE.at(i).second;
      }
    }

    for(int i = 0; i < trackEIn.size(); ++i){
      int eParent = trackEIn.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if parent = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0InDep += trackEIn.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMInDep += trackEIn.at(i).second;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMInDep += trackEIn.at(i).second;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMInDep += trackEIn.at(i).second;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMInDep += trackE.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMInDep += trackE.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMInDep += trackE.at(i).second;
      }

    }

    for(int i = 0; i < trackEOut.size(); ++i){
      int eParent = trackEOut.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if paremt = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0OutDep += trackEOut.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMOutDep += trackEOut.at(i).second;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMOutDep += trackEOut.at(i).second;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMOutDep += trackEOut.at(i).second;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMOutDep += trackE.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMOutDep += trackE.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMOutDep += trackE.at(i).second;
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
}*/

void DunePrismAnalyzer::AnalyzeStops(){
  for (int ie = 0; ie < nEntries; ++ie){
    if(!(ie%100))std::cout<<"Entry " << ie<<std::endl;
    inTree->GetEntry(ie);
    if(nuPID != 14) continue;
//    #ifdef DEBUG
//    std::cout << i << std::endl;
//    std::cout<<Enu<<std::endl;
//    std::cout<<nuPID<<std::endl;
//    #endif 
 
    //To Do: make a vector of trees
    //that correspond to different 
    //stops off axis
    
    //////////////Vertex in Detector?
    #ifdef DEBUG
//    std::cout << "vtxX " << EvtVtx[0] << std::endl;
//    std::cout << "vtxY " << fabs(EvtVtx[1]) << std::endl; 
//    std::cout << "vtxZ " << fabs(EvtVtx[2]) << std::endl; 
    #endif 

    int stop_num;
    bool noStops = true;

    double wallX[2] = {0.,0.};
    double wallY = 0.;
    double wallZ = 0.;

    for(int i = 0; i < nStops; ++i){

      
      wallX[0] = (shift.at(i) - dimension.at(i)[0]*.5 + FV.at(i)[0]);
//      std::cout <<"WALLLL " << wallX[0] << std::endl;
      wallX[1] = (shift.at(i) + dimension.at(i)[0]*.5 - FV.at(i)[0]);
      wallY = dimension.at(i)[1]*.5 - FV.at(i)[1];
      wallZ = dimension.at(i)[2]*.5 - FV.at(i)[2];

      //First: check if in y & z FV
      if( fabs(EvtVtx[1]) > wallY ) continue;
      else if( fabs(EvtVtx[2]) > wallZ ) continue;

      //Check if in x 
      else if( EvtVtx[0] > wallX[1] || EvtVtx[0] < wallX[0]) continue;

      stop_num = i;
      noStops = false;

      #ifdef DEBUG
//      std::cout << "Inside FV "<< stop_num << std::endl;
      #endif

      break;//Leave after finding first stop
            //that contains the event
    }

    if(noStops) continue;
    /////////////////////////////////////
    
    Enu.at(stop_num) = ekina;
    vtx_X.at(stop_num) = EvtVtx[0];
    vtx_Y.at(stop_num) = EvtVtx[1];
    vtx_Z.at(stop_num) = EvtVtx[2];
    eventNum.at(stop_num) = ev; 
   
    bool found = false;
    bool muEnd = false;

    flagExitBack.at(stop_num) = 0;
    flagExitFront.at(stop_num) = 0;
    flagExitY.at(stop_num) = 0;
    flagExitXHigh.at(stop_num) = 0;
    flagExitXLow.at(stop_num) = 0;
    flagMuContained.at(stop_num) = 0;
    flagNoEHadOut.at(stop_num) = 0;

    muExitingPX.at(stop_num) = 0.;
    muExitingPY.at(stop_num) = 0.;
    muExitingPZ.at(stop_num) = 0.;

    eHadOutDep.at(stop_num) = 0.;
    eHadInDep.at(stop_num) = 0.;
    eHadTotalDep.at(stop_num) = 0.;

    ePiCTotalDep.at(stop_num) = 0.; 
    ePiCInDep.at(stop_num) = 0.;
    ePiCOutDep.at(stop_num) = 0.;

    eProtonTotalDep.at(stop_num) = 0.; 
    eProtonInDep.at(stop_num) = 0.;
    eProtonOutDep.at(stop_num) = 0.;

    eMuDep.at(stop_num) = 0.;
    eMuTotalDep.at(stop_num) = 0.;

    std::vector<int>::iterator itPi0;
    std::map<int,int>::iterator itGamma;
    std::vector<int>::iterator itPiC;
    std::vector<int>::iterator itProton;
    std::vector<int>::iterator itMu;

    ///New
    trackMu.clear();
    trackPiC.clear();
    trackProton.clear();
    
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
          itMu = std::find( trackMu.begin(), trackMu.end(), track[i] );
          if(itMu == trackMu.end() ){
            trackMu.push_back(track[i]);
          }
        }                
        
        eMuDep.at(stop_num) += edep[i];
        eMuTotalDep.at(stop_num) += edep[i];

        //check if it hits the walls
        if(xe[i] < wallX[0] - FV.at(stop_num)[0]){
          flagExitXLow.at(stop_num) = true; 
          muExitingPX.at(stop_num) = pxe[i];
          muExitingPY.at(stop_num) = pye[i];
          muExitingPZ.at(stop_num) = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted XLow "<< 
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
        else if(xe[i] > wallX[1] + FV.at(stop_num)[0]){ 
          flagExitXHigh.at(stop_num) = true; 
          muExitingPX.at(stop_num) = pxe[i];
          muExitingPY.at(stop_num) = pye[i];
          muExitingPZ.at(stop_num) = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted right "<< 
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }

        if(abs(ye[i]) > wallY + FV.at(stop_num)[1]){ 
          flagExitY.at(stop_num) = true; 
          muExitingPX.at(stop_num) = pxe[i];
          muExitingPY.at(stop_num) = pye[i];
          muExitingPZ.at(stop_num) = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted Y "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }

        if(ze[i] > wallZ + FV.at(stop_num)[2]){ 
          flagExitBack.at(stop_num) = true; 
          muExitingPX.at(stop_num) = pxe[i];
          muExitingPY.at(stop_num) = pye[i];
          muExitingPZ.at(stop_num) = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted back "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
        else if(ze[i] < -wallZ - FV.at(stop_num)[2]){ 
          flagExitFront.at(stop_num) = true; 
          muExitingPX.at(stop_num) = pxe[i];
          muExitingPY.at(stop_num) = pye[i];
          muExitingPZ.at(stop_num) = pze[i];
          #ifdef MUEXITDEBUG
          std::cout<<"Muon exitted front "<<
          xe[i] << " "<<ye[i] << " "<<ze[i]<<
          std::endl;
          #endif
          muEnd = true;
        }
      }

      else if(PID[i] == 13 && muEnd){//Energy dep outside detector
        eMuTotalDep.at(stop_num) += edep[i];
      }

      //Look for hadronic energy deposits
      else if(PID[i] == 2212){
        itProton = std::find( trackProton.begin(), trackProton.end(), track[i] );
        if(itProton == trackProton.end() ){
          trackProton.push_back(track[i]);
        }       
        eProtonTotalDep.at(stop_num) += edep[i];

        //Inside walls. Excuse terrible naming
        if(xe[i] > (wallX[0] - FV.at(stop_num)[0]) && xe[i] < (wallX[1] + FV.at(stop_num)[0])&&
          abs(ye[i]) < (wallY + FV.at(stop_num)[0]) && abs(ze[i]) < (wallZ + FV.at(stop_num)[0])
          ){
        
          //Between wall and FV
          if(xe[i] < wallX[0] || xe[i] > wallX[1] ||
            abs(ye[i]) > wallY || abs(ze[i]) > wallZ){
            
            eProtonOutDep.at(stop_num) += edep[i]; 
          }

          else eProtonInDep.at(stop_num) += edep[i];
        }     
      }
      else if(abs(PID[i]) == 211){
        itPiC = std::find( trackPiC.begin(), trackPiC.end(), track[i] );
        if(itPiC == trackPiC.end() ){
          trackPiC.push_back(track[i]);
        }
        ePiCTotalDep.at(stop_num) += edep[i];

        //Inside walls. Excuse terrible naming
        if(xe[i] > (wallX[0] - FV.at(stop_num)[0]) && xe[i] < (wallX[1] + FV.at(stop_num)[0])&&
          abs(ye[i]) < (wallY + FV.at(stop_num)[0]) && abs(ze[i]) < (wallZ + FV.at(stop_num)[0])
          ){
        
          //Between wall and FV
          if(xe[i] < wallX[0] || xe[i] > wallX[1] ||
            abs(ye[i]) > wallY || abs(ze[i]) > wallZ){
            
            ePiCOutDep.at(stop_num) += edep[i]; 
          }

          else ePiCInDep.at(stop_num) += edep[i];
        }
      }

      //Look for Pi0. Add in track id if not in vector
      else if(PID[i] == 111){
        itPi0 = std::find( trackPi0.begin(), trackPi0.end(), track[i] );
        if(itPi0 == trackPi0.end() ){
          trackPi0.push_back(track[i]);
        }

      }

      //Look for gamma. Add in track and parid if not in map
      else if(PID[i] == 22){
        if(trackGamma.find(track[i]) == trackGamma.end() ){
          trackGamma[track[i]] = parid[i];
        }
      }

      //Add in all edeps and parids from electrons
      else if(abs(PID[i]) == 11){        
        trackE.push_back( std::make_pair(parid[i],edep[i]) );

        //Inside walls. Excuse terrible naming
        if(xe[i] > (wallX[0] - FV.at(stop_num)[0]) && xe[i] < (wallX[1] + FV.at(stop_num)[0])&&
          abs(ye[i]) < (wallY + FV.at(stop_num)[0]) && abs(ze[i]) < (wallZ + FV.at(stop_num)[0])
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

    flagMuContained.at(stop_num) = !(flagExitFront.at(stop_num) || flagExitBack.at(stop_num) ||
                         flagExitXLow.at(stop_num) || flagExitXHigh.at(stop_num) ||
                         flagExitY.at(stop_num));
        
    if(trackE.size() == 0)std::cout<<"found no e+-"<<std::endl;
    
    ePi0TotalDep.at(stop_num) = 0.; 
    ePi0InDep.at(stop_num) = 0.;
    ePi0OutDep.at(stop_num) = 0.;

    ePiCEMTotalDep.at(stop_num) = 0.; 
    ePiCEMInDep.at(stop_num) = 0.;
    ePiCEMOutDep.at(stop_num) = 0.;

    eProtonEMTotalDep.at(stop_num) = 0.; 
    eProtonEMInDep.at(stop_num) = 0.;
    eProtonEMOutDep.at(stop_num) = 0.;

    eMuEMTotalDep.at(stop_num) = 0.; 
    eMuEMInDep.at(stop_num) = 0.;
    eMuEMOutDep.at(stop_num) = 0.;

    eResidualEMTotalDep.at(stop_num) = 0.;
    eResidualEMInDep.at(stop_num) = 0.;
    eResidualEMOutDep.at(stop_num) = 0.;

    //Iterate over e+- edeps, check if parid in gamma map
    //if so: check if gamma parid in pi0 vector
    for(int i = 0; i < trackE.size(); ++i){
      int eParent = trackE.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if parent = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0TotalDep.at(stop_num) += trackE.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMTotalDep.at(stop_num) += trackE.at(i).second;
//          std::cout << "Found pic->gamma->e"<<std::endl;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMTotalDep.at(stop_num) += trackE.at(i).second;
//          std::cout << "Found proton->gamma->e"<<std::endl;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMTotalDep.at(stop_num) += trackE.at(i).second;
 //         std::cout << "Found mu->gamma->e"<<std::endl;
        }
        else{
          eResidualEMTotalDep.at(stop_num) += trackE.at(i).second;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMTotalDep.at(stop_num) += trackE.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMTotalDep.at(stop_num) += trackE.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMTotalDep.at(stop_num) += trackE.at(i).second;
      }
      else{
        eResidualEMTotalDep.at(stop_num) += trackE.at(i).second;
      }
    }

    for(int i = 0; i < trackEIn.size(); ++i){
      int eParent = trackEIn.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if parent = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0InDep.at(stop_num) += trackEIn.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMInDep.at(stop_num) += trackEIn.at(i).second;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMInDep.at(stop_num) += trackEIn.at(i).second;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMInDep.at(stop_num) += trackEIn.at(i).second;
        }
        else{
          eResidualEMInDep.at(stop_num) += trackEIn.at(i).second;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMInDep.at(stop_num) += trackEIn.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMInDep.at(stop_num) += trackEIn.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMInDep.at(stop_num) += trackEIn.at(i).second;
      }
      else{
        eResidualEMInDep.at(stop_num) += trackEIn.at(i).second;
      }

    }

    for(int i = 0; i < trackEOut.size(); ++i){
      int eParent = trackEOut.at(i).first;

      itPiC = std::find(trackPiC.begin(), trackPiC.end(),eParent);
      itProton = std::find(trackProton.begin(), trackProton.end(),eParent);
      itMu = std::find(trackMu.begin(), trackMu.end(),eParent);

      // if paremt = gamma
      if(trackGamma.find(eParent) != trackGamma.end()){
        int gammaParent = trackGamma[eParent];
        itPi0 = std::find(trackPi0.begin(),trackPi0.end(),gammaParent);
        if(itPi0 != trackPi0.end()){//Successfully found pi0->g->e 
          ePi0OutDep.at(stop_num) += trackEOut.at(i).second; //add in edep
        }

        else if(std::find(trackPiC.begin(),trackPiC.end(),gammaParent)
                != trackPiC.end()){
          ePiCEMOutDep.at(stop_num) += trackEOut.at(i).second;
        }

        else if(std::find(trackProton.begin(),trackProton.end(),gammaParent)
                != trackProton.end()){
          eProtonEMOutDep.at(stop_num) += trackEOut.at(i).second;
        }
        
        else if(std::find(trackMu.begin(),trackMu.end(),gammaParent)
                != trackMu.end()){
          eMuEMOutDep.at(stop_num) += trackEOut.at(i).second;
        }
        else{
          eResidualEMOutDep.at(stop_num) += trackEOut.at(i).second;
        }
      }

      // if parent = charged pion
      else if(itPiC != trackPiC.end()){
        ePiCEMOutDep.at(stop_num) += trackEOut.at(i).second;
      }

      // if parent = proton
      else if(itProton != trackProton.end()){
        eProtonEMOutDep.at(stop_num) += trackEOut.at(i).second;
      }

      // if parent = muon
      else if(itMu != trackMu.end()){
        eMuEMOutDep.at(stop_num) += trackEOut.at(i).second;
      }
      else{
        eResidualEMOutDep.at(stop_num) += trackEOut.at(i).second;
      }

    }

    flagNoEHadOut.at(stop_num) = (ePi0OutDep.at(stop_num) == 0 && ( ePi0TotalDep.at(stop_num) == ePi0InDep.at(stop_num) )
                     && eHadOutDep.at(stop_num) == 0 
                     && ( eHadTotalDep.at(stop_num) == eHadInDep.at(stop_num) ) );
    ////////////////////////////////////

    eHadTrueCharged.at(stop_num) = 0.;
    eHadTrueTotal.at(stop_num) = 0.;
    eMuTrue.at(stop_num) = 0.;

    nMu.at(stop_num) = 0;
    nPi0.at(stop_num) = 0;
    nPiC.at(stop_num) = 0;
    nProton.at(stop_num) = 0;
    nNeutron.at(stop_num) = 0;

    //Start of particle loop
//    std::cout << "Starting true FS particle loop "<< ni << std::endl;

    for(int ip = 0; ip < ni; ++ip){
      //std::cout << ip << std::endl;
      if(PIDi[ip] == 13){
        eMuTrue.at(stop_num) += ekini[ip] + mi[ip];
        nMu.at(stop_num)++;
      }
      else if(PIDi[ip] == 111){
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        nPi0.at(stop_num)++;
      }
      else if(abs(PIDi[ip]) == 211 ){
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        nPiC.at(stop_num)++;
      }
      else if(PIDi[ip] == 2212){
        eHadTrueCharged.at(stop_num) += ekini[ip];
        nProton.at(stop_num)++;
      }
      else if(PIDi[ip] == 2112){
        eHadTrueTotal.at(stop_num) += ekini[ip];
        nNeutron.at(stop_num)++;
      }
    }
    eHadTrueTotal.at(stop_num) += eHadTrueCharged.at(stop_num); 
    ///////////////////////////////////

    stopsTrees.at(stop_num)->Fill();
  }
}

/*void DunePrismAnalyzer::SetBranches(){

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

  eventTree->Branch("ePiCEMTotalDep",&ePiCEMTotalDep,"ePiCEMTotalDep/D");
  eventTree->Branch("ePiCEMInDep",&ePiCEMInDep,"ePiCEMInDep/D");
  eventTree->Branch("ePiCEMOutDep",&ePiCEMOutDep,"ePiCEMOutDep/D");

  eventTree->Branch("eProtonEMTotalDep",&eProtonEMTotalDep,"eProtonEMTotalDep/D");
  eventTree->Branch("eProtonEMInDep",&eProtonEMInDep,"eProtonEMInDep/D");
  eventTree->Branch("eProtonEMOutDep",&eProtonEMOutDep,"eProtonEMOutDep/D");

  eventTree->Branch("eMuEMTotalDep",&eMuEMTotalDep,"eMuEMTotalDep/D");
  eventTree->Branch("eMuEMInDep",&eMuEMInDep,"eMuEMInDep/D");
  eventTree->Branch("eMuEMOutDep",&eMuEMOutDep,"eMuEMOutDep/D");

  eventTree->Branch("eHadTrueCharged",&eHadTrueCharged,"eHadTrueCharged/D");
  eventTree->Branch("eHadTrueTotal",&eHadTrueTotal,"eHadTrueTotal/D");
  eventTree->Branch("eMuTrue",&eMuTrue,"eMuTrue/D");
  eventTree->Branch("nMu",&nMu,"nMu/I");
  eventTree->Branch("nPi0",&nPi0,"nPi0/I");
  eventTree->Branch("nPiC",&nPiC,"nPiC/I");
  eventTree->Branch("nProton",&nProton,"nProton/I");
  eventTree->Branch("nNeutron",&nNeutron,"nNeutron/I");

}*/

void DunePrismAnalyzer::SetInBranches(){

  inTree->SetBranchAddress("ni",&ni);
  inTree->SetBranchAddress("pidi",&PIDi);
  inTree->SetBranchAddress("pxi",&pxi);
  inTree->SetBranchAddress("pyi",&pyi);
  inTree->SetBranchAddress("pzi",&pzi);
  inTree->SetBranchAddress("ekini",&ekini);
  inTree->SetBranchAddress("mi",&mi);

  inTree->SetBranchAddress("ev",&ev);
  inTree->SetBranchAddress("ekina",&ekina);
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
}

void DunePrismAnalyzer::SetOutBranches(int stop){
//  InitVarsStop(stop);
  stopsTrees.at(stop)->Branch("eventNum",&eventNum.at(stop));
  stopsTrees.at(stop)->Branch("Enu",&Enu.at(stop),"Enu/D");
  stopsTrees.at(stop)->Branch("vtx_X",&vtx_X.at(stop),"vtx_X/D");
  stopsTrees.at(stop)->Branch("vtx_Y",&vtx_Y.at(stop),"vtx_Y/D");
  stopsTrees.at(stop)->Branch("vtx_Z",&vtx_Z.at(stop),"vtx_Z/D");

  stopsTrees.at(stop)->Branch("flagExitBack",&flagExitBack.at(stop));
  stopsTrees.at(stop)->Branch("flagExitFront",&flagExitFront.at(stop));
  stopsTrees.at(stop)->Branch("flagExitY",&flagExitY.at(stop));
  stopsTrees.at(stop)->Branch("flagExitXLow",&flagExitXLow.at(stop));
  stopsTrees.at(stop)->Branch("flagExitXHigh",&flagExitXHigh.at(stop));
  stopsTrees.at(stop)->Branch("flagNoEHadOut",&flagNoEHadOut.at(stop));
  stopsTrees.at(stop)->Branch("flagMuContained",&flagMuContained.at(stop));

  stopsTrees.at(stop)->Branch("muExitingPX",&muExitingPX.at(stop),"muExitingPX/D");
  stopsTrees.at(stop)->Branch("muExitingPY",&muExitingPY.at(stop),"muExitingPY/D");
  stopsTrees.at(stop)->Branch("muExitingPZ",&muExitingPZ.at(stop),"muExitingPZ/D");

  stopsTrees.at(stop)->Branch("eHadOutDep",&eHadOutDep.at(stop),"eHadOutDep/D");
  stopsTrees.at(stop)->Branch("eHadInDep",&eHadInDep.at(stop),"eHadInDep/D");
  stopsTrees.at(stop)->Branch("eHadTotalDep",&eHadTotalDep.at(stop),"eHadTotalDep/D");

  stopsTrees.at(stop)->Branch("eMuDep",&eMuDep.at(stop),"eMuDep/D");
  stopsTrees.at(stop)->Branch("eMuTotalDep",&eMuTotalDep.at(stop),"eMuTotalDep/D");

  stopsTrees.at(stop)->Branch("ePi0TotalDep",&ePi0TotalDep.at(stop),"ePi0TotalDep/D");
  stopsTrees.at(stop)->Branch("ePi0InDep",&ePi0InDep.at(stop),"ePi0InDep/D");
  stopsTrees.at(stop)->Branch("ePi0OutDep",&ePi0OutDep.at(stop),"ePi0OutDep/D");

  stopsTrees.at(stop)->Branch("ePiCTotalDep",&ePiCTotalDep.at(stop),"ePiCTotalDep/D");
  stopsTrees.at(stop)->Branch("ePiCInDep",&ePiCInDep.at(stop),"ePiCInDep/D");
  stopsTrees.at(stop)->Branch("ePiCOutDep",&ePiCOutDep.at(stop),"ePiCOutDep/D");

  stopsTrees.at(stop)->Branch("eProtonTotalDep",&eProtonTotalDep.at(stop),"eProtonTotalDep/D");
  stopsTrees.at(stop)->Branch("eProtonInDep",&eProtonInDep.at(stop),"eProtonInDep/D");
  stopsTrees.at(stop)->Branch("eProtonOutDep",&eProtonOutDep.at(stop),"eProtonOutDep/D");

  stopsTrees.at(stop)->Branch("ePiCEMTotalDep",&ePiCEMTotalDep.at(stop),"ePiCEMTotalDep/D");
  stopsTrees.at(stop)->Branch("ePiCEMInDep",&ePiCEMInDep.at(stop),"ePiCEMInDep/D");
  stopsTrees.at(stop)->Branch("ePiCEMOutDep",&ePiCEMOutDep.at(stop),"ePiCEMOutDep/D");

  stopsTrees.at(stop)->Branch("eProtonEMTotalDep",&eProtonEMTotalDep.at(stop),"eProtonEMTotalDep/D");
  stopsTrees.at(stop)->Branch("eProtonEMInDep",&eProtonEMInDep.at(stop),"eProtonEMInDep/D");
  stopsTrees.at(stop)->Branch("eProtonEMOutDep",&eProtonEMOutDep.at(stop),"eProtonEMOutDep/D");

  stopsTrees.at(stop)->Branch("eMuEMTotalDep",&eMuEMTotalDep.at(stop),"eMuEMTotalDep/D");
  stopsTrees.at(stop)->Branch("eMuEMInDep",&eMuEMInDep.at(stop),"eMuEMInDep/D");
  stopsTrees.at(stop)->Branch("eMuEMOutDep",&eMuEMOutDep.at(stop),"eMuEMOutDep/D");

  stopsTrees.at(stop)->Branch("eResidualEMTotalDep",&eResidualEMTotalDep.at(stop),"eResidualEMTotalDep/D");
  stopsTrees.at(stop)->Branch("eResidualEMInDep",&eResidualEMInDep.at(stop),"eResidualEMInDep/D");
  stopsTrees.at(stop)->Branch("eResidualEMOutDep",&eResidualEMOutDep.at(stop),"eResidualEMOutDep/D");

  stopsTrees.at(stop)->Branch("eHadTrueCharged",&eHadTrueCharged.at(stop),"eHadTrueCharged/D");
  stopsTrees.at(stop)->Branch("eHadTrueTotal",&eHadTrueTotal.at(stop),"eHadTrueTotal/D");
  stopsTrees.at(stop)->Branch("eMuTrue",&eMuTrue.at(stop),"eMuTrue/D");
  stopsTrees.at(stop)->Branch("nMu",&nMu.at(stop),"nMu/I");
  stopsTrees.at(stop)->Branch("nPi0",&nPi0.at(stop),"nPi0/I");
  stopsTrees.at(stop)->Branch("nPiC",&nPiC.at(stop),"nPiC/I");
  stopsTrees.at(stop)->Branch("nProton",&nProton.at(stop),"nProton/I");
  stopsTrees.at(stop)->Branch("nNeutron",&nNeutron.at(stop),"nNeutron/I");

}

void DunePrismAnalyzer::InitVarsStop(){
  eventNum.push_back(0);
  Enu.push_back(0.);
//  std::array<double,3> vtx = {0.,0.,0.};
  vtx_X.push_back(0.);
  vtx_Y.push_back(0.);
  vtx_Z.push_back(0.);
  flagExitBack.push_back(0);
  flagExitFront.push_back(0);
  flagExitY.push_back(0);
  flagExitXHigh.push_back(0);
  flagExitXLow.push_back(0);
  flagNoEHadOut.push_back(0);
  flagMuContained.push_back(0);

  muExitingPX.push_back(0.);
  muExitingPY.push_back(0.);
  muExitingPZ.push_back(0.);

  eHadOutDep.push_back(0.);
  eHadInDep.push_back(0.);
  eHadTotalDep.push_back(0.);

  ePi0TotalDep.push_back(0.);
  ePi0InDep.push_back(0.);
  ePi0OutDep.push_back(0.);

  ePiCTotalDep.push_back(0.);
  ePiCInDep.push_back(0.);
  ePiCOutDep.push_back(0.);

  eProtonTotalDep.push_back(0.);
  eProtonInDep.push_back(0.);
  eProtonOutDep.push_back(0.);

  ePiCEMTotalDep.push_back(0.);
  ePiCEMInDep.push_back(0.);
  ePiCEMOutDep.push_back(0.);

  eProtonEMTotalDep.push_back(0.);
  eProtonEMInDep.push_back(0.);
  eProtonEMOutDep.push_back(0.);

  eMuEMTotalDep.push_back(0.);
  eMuEMInDep.push_back(0.);
  eMuEMOutDep.push_back(0.);

  eResidualEMTotalDep.push_back(0.);
  eResidualEMInDep.push_back(0.);
  eResidualEMOutDep.push_back(0.);

  eMuDep.push_back(0.);
  eMuTotalDep.push_back(0.);

  eHadTrueCharged.push_back(0.);
  eHadTrueTotal.push_back(0.);
  eMuTrue.push_back(0.);

  nMu.push_back(0.);
  nPi0.push_back(0.);
  nPiC.push_back(0.);
  nProton.push_back(0.);
  nNeutron.push_back(0.);

}
void DunePrismAnalyzer::Finalize(){
  fout->cd();
  eventTree->Write();
  //gOutTree->Write();
  
  fout->Close();
  fin->Close();

}

void DunePrismAnalyzer::FinalizeStops(){
  fout->cd();
  for(int i = 0; i < stopsTrees.size(); ++i){
    stopsTrees.at(i)->Write();  
  }
  
  fout->Close();
  fin->Close();

}
