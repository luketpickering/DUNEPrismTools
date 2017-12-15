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

void DunePrismAnalyzer::AnalyzeStops(){
  for (int ie = 0; ie < nEntries; ++ie){
    if(!(ie%100))std::cout<<"Entry " << ie<<std::endl;
    inTree->GetEntry(ie);
    if(nuPID != 14) continue;

    int stop_num;
    bool noStops = true;

    double wallX[2] = {0.,0.};
    double wallY = 0.;
    double wallZ = 0.;
 

    for(int i = 0; i < nStops; ++i){

      
      wallX[0] = (shift.at(i) - dimension.at(i)[0]*.5 + FV.at(i)[0]);
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

    double muonBoundX[2] = {
      wallX[0] - FV.at(stop_num)[0],
      wallX[1] + FV.at(stop_num)[0]
    };
    double muonBoundY = wallY + FV.at(stop_num)[1];
    double muonBoundZ = wallZ + FV.at(stop_num)[1];
   
    bool found = false;
    bool muEnd = false;

    eHadTrueCharged.at(stop_num) = 0.;
    eHadTrueTotal.at(stop_num) = 0.;
    eMuTrue.at(stop_num) = 0.;
    eGammaTrue.at(stop_num) = 0.;

    pMuTrueX.at(stop_num) = 0.;
    pMuTrueY.at(stop_num) = 0.;
    pMuTrueZ.at(stop_num) = 0.;

    nMu.at(stop_num) = 0;
    nPi0.at(stop_num) = 0;
    nPiC.at(stop_num) = 0;
    nProton.at(stop_num) = 0;
    nNeutron.at(stop_num) = 0;
    nGamma.at(stop_num) = 0;

    std::map<int,int> chain;

    int nBindino = 0;
    int nOther = 0;

    FSHadrons.clear();

    //Start of particle loop
    for(int ip = 0; ip < ni; ++ip){
      //std::cout << ip << std::endl;
      if (PIDi[ip] == 2000000101){
        nBindino++;
        continue;//Skip bindino
      }
      chain[ip + 1] = 0; 
      if(PIDi[ip] == 13){ //Initial muon 
 //       std::cout << "Found " << PIDi[ip] << " " << ip + 1 << std::endl; 
        FSMuon = new DepoMuon(
          PIDi[ip], (ip + 1), ekini[ip],
          muonBoundX, muonBoundY, muonBoundZ); 

        eMuTrue.at(stop_num) = ekini[ip] + mi[ip];
        pMuTrueX.at(stop_num) = pxi[ip];
        pMuTrueY.at(stop_num) = pyi[ip];
        pMuTrueZ.at(stop_num) = pzi[ip];

        Q2True.at(stop_num) = fabs( pow( (ekina - eMuTrue.at(stop_num)), 2 ) 
                                 - pow( (pxi[ip] - pxa), 2 )
                                 - pow( (pyi[ip] - pya), 2 )
                                 - pow( (pzi[ip] - pza), 2 ) );

        yTrue.at(stop_num) = 1  - ekini[ip]/ekina;
        double mN = .93827208;
        W_rest.at(stop_num) = sqrt(-Q2True.at(stop_num) + 2 * mN * (ekina - eMuTrue.at(stop_num)) + mN * mN);

        nMu.at(stop_num)++;
      }
      else if(PIDi[ip] == 111){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        nPi0.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(abs(PIDi[ip]) == 211 ){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        nPiC.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(PIDi[ip] == 2212){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip];
        nProton.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(PIDi[ip] == 2112){
 //       std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueTotal.at(stop_num) += ekini[ip];
        nNeutron.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(PIDi[ip] == 22){
  //      std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eGammaTrue.at(stop_num) += ekini[ip];
        nGamma.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
        //Add it to hadrons because I'm lazy 
        //to make a new thing
      }
      else{
//        std::cout << "Found Other! " << PIDi[ip] <<  " " << ip + 1 << std::endl;
        nOther++;
      }
    }

    int nPrimary = ni - nBindino - nOther;
//    std::cout << "NHadrons: " << FSHadrons.size() << std::endl;
    if (nMu.at(stop_num) == 0) continue;

    eHadTrueTotal.at(stop_num) += eHadTrueCharged.at(stop_num); 
    ///////////////////////////////////

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

    eMuPrimaryDep.at(stop_num) = 0.;
    eMuSecondaryDep.at(stop_num) = 0.;

    eHadPrimaryDepIn.at(stop_num) = 0.;
    eHadSecondaryDepIn.at(stop_num) = 0.;
    eHadPrimaryDepOut.at(stop_num) = 0.;
    eHadSecondaryDepOut.at(stop_num) = 0.;

    eProtonPrimaryDepIn.at(stop_num) = 0.;
    eProtonSecondaryDepIn.at(stop_num) = 0.;
    eProtonPrimaryDepOut.at(stop_num) = 0.;
    eProtonSecondaryDepOut.at(stop_num) = 0.;

    eNeutronPrimaryDepIn.at(stop_num) = 0.;
    eNeutronSecondaryDepIn.at(stop_num) = 0.;
    eNeutronPrimaryDepOut.at(stop_num) = 0.;
    eNeutronSecondaryDepOut.at(stop_num) = 0.;

    ePiCPrimaryDepIn.at(stop_num) = 0.;
    ePiCSecondaryDepIn.at(stop_num) = 0.;
    ePiCPrimaryDepOut.at(stop_num) = 0.;
    ePiCSecondaryDepOut.at(stop_num) = 0.;

    eGammaPrimaryDepIn.at(stop_num) = 0.;
    eGammaSecondaryDepIn.at(stop_num) = 0.;
    eGammaPrimaryDepOut.at(stop_num) = 0.;
    eGammaSecondaryDepOut.at(stop_num) = 0.;

    ePi0PrimaryDepIn.at(stop_num) = 0.;
    ePi0SecondaryDepIn.at(stop_num) = 0.;
    ePi0PrimaryDepOut.at(stop_num) = 0.;
    ePi0SecondaryDepOut.at(stop_num) = 0.;

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
//    std::cout << "Starting steps " << std::endl;
    for(int i = 0; i < nstep; ++i){
      if(chain.find(track[i]) == chain.end()){//Not in chain
      
        if(chain.find(parid[i]) == chain.end()){//Error
          std::cout << "ERROR: PARENT NOT FOUND" << std::endl;
          break;
        }
        chain[track[i]] = parid[i];
      }
    
      if(chain[track[i]] == 0){//Primary

        if(PID[i] == 13){//Muon
          if(FSMuon->trackID != track[i]){
            std::cout << "ERROR: Wrong muon track" << std::endl;  
          }

          FSMuon->eDepPrimary += edep[i];
          FSMuon->xf = xe[i];
          FSMuon->yf = ye[i];
          FSMuon->zf = ze[i];

          FSMuon->pxf = pxe[i];
          FSMuon->pyf = pye[i];
          FSMuon->pzf = pze[i];
        }
        else{//Hadron
          if(!(PID[i] == 2212 || PID[i] == 2112 || abs(PID[i]) == 211 ||
               PID[i] == 111  || PID[i] ==   22))continue;
               
          if(xe[i] <= FSHadrons[track[i]]->xBound[0] ||
             xe[i] >= FSHadrons[track[i]]->xBound[1] ||
              fabs(ye[i]) >= FSHadrons[track[i]]->yBound ||
              fabs(ze[i]) >= FSHadrons[track[i]]->zBound){
            
            FSHadrons[track[i]]->eDepPrimaryOut += edep[i];
          }
          else{
            FSHadrons[track[i]]->eDepPrimaryIn += edep[i];
          }         
        }
      }
      else{//Secondary

        //Find ultimate parent
        int itChain = track[i];
        while (chain[itChain] != 0){
          itChain = chain[itChain];                
        }
        if(itChain == 1){//Muon
          FSMuon->eDepSecondary += edep[i];
        }
        else{//Hadron
          if(!(PID[i] == 2212 || PID[i] == 2112 || abs(PID[i]) == 211 ||
               PID[i] == 111  || PID[i] ==   22))continue;
          if(FSHadrons.find(itChain) == FSHadrons.end()){
//            std::cout << "ERROR: ULTIMATE HADRON NOT IN FS " << itChain << std::endl;
            continue;
          }
          if(xe[i] <= FSHadrons[itChain]->xBound[0] ||
             xe[i] >= FSHadrons[itChain]->xBound[1] ||
              fabs(ye[i]) >= FSHadrons[itChain]->yBound ||
              fabs(ze[i]) >= FSHadrons[itChain]->zBound){
            
            FSHadrons[itChain]->eDepSecondaryOut += edep[i];
          }
          else{
            FSHadrons[itChain]->eDepSecondaryIn += edep[i];
          }
        }
      }
    }
//    std::cout << "Ending steps " << std::endl;
    eMuPrimaryDep.at(stop_num) += FSMuon->eDepPrimary;
    eMuSecondaryDep.at(stop_num) += FSMuon->eDepSecondary;

    FSMuon->CheckContained();

    flagExitXLow.at(stop_num) = FSMuon->flagExitXLow;
    flagExitXHigh.at(stop_num) = FSMuon->flagExitXHigh;
    flagExitY.at(stop_num) = FSMuon->flagExitY;
    flagExitBack.at(stop_num) = FSMuon->flagExitBack;
    flagExitFront.at(stop_num) = FSMuon->flagExitFront;
    flagMuContained.at(stop_num) = FSMuon->flagMuContained;

    muExitingPX.at(stop_num) = FSMuon->pxf;
    muExitingPY.at(stop_num) = FSMuon->pyf;
    muExitingPZ.at(stop_num) = FSMuon->pzf;   

    flagNoEHadOut.at(stop_num) = 0;
    
    std::map<int,DepoHadron*>::iterator itHad;
    for(itHad = FSHadrons.begin(); itHad != FSHadrons.end(); ++itHad){
      eHadPrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
      eHadSecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

      eHadPrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
      eHadSecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;
       
      //Splitting energy up into various particles
      switch( abs(itHad->second->PDG) ){
        case 2212: 
          eProtonPrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
          eProtonSecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

          eProtonPrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
          eProtonSecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;

          break;
        case 2112:
          eNeutronPrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
          eNeutronSecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

          eNeutronPrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
          eNeutronSecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;

          break;
        case  211:
          ePiCPrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
          ePiCSecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

          ePiCPrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
          ePiCSecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;

          break;
        case   22:
          eGammaPrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
          eGammaSecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

          eGammaPrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
          eGammaSecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;

          break;
        case  111:
          ePi0PrimaryDepIn.at(stop_num) += itHad->second->eDepPrimaryIn;
          ePi0SecondaryDepIn.at(stop_num) += itHad->second->eDepSecondaryIn;

          ePi0PrimaryDepOut.at(stop_num) += itHad->second->eDepPrimaryOut;
          ePi0SecondaryDepOut.at(stop_num) += itHad->second->eDepSecondaryOut;

          break;
      }
      //

      if(itHad->second->eDepPrimaryOut + itHad->second->eDepSecondaryOut > 0.000000000000000001){
        flagNoEHadOut.at(stop_num) = 1;
      }

    }

    eReco.at(stop_num) = eMuPrimaryDep.at(stop_num) + eMuSecondaryDep.at(stop_num) + eHadPrimaryDepIn.at(stop_num) + eHadSecondaryDepIn.at(stop_num);


//    std::cout << eReco.at(stop_num) << std::endl;

    ///////////////////////////////
    stopsTrees.at(stop_num)->Fill();
  }
}


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
  inTree->SetBranchAddress("pxa",&pxa);
  inTree->SetBranchAddress("pya",&pya);
  inTree->SetBranchAddress("pza",&pza);
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

  stopsTrees.at(stop)->Branch("eMuPrimaryDep",&eMuPrimaryDep.at(stop),"eMuPrimaryDep/D");
  stopsTrees.at(stop)->Branch("eMuSecondaryDep",&eMuSecondaryDep.at(stop),"eMuSecondaryDep/D");

  stopsTrees.at(stop)->Branch("eHadPrimaryDepIn",&eHadPrimaryDepIn.at(stop),"eHadPrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("eHadSecondaryDepIn",&eHadSecondaryDepIn.at(stop),"eHadSecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("eHadPrimaryDepOut",&eHadPrimaryDepOut.at(stop),"eHadPrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("eHadSecondaryDepOut",&eHadSecondaryDepOut.at(stop),"eHadSecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("eProtonPrimaryDepIn",&eProtonPrimaryDepIn.at(stop),"eProtonPrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("eProtonSecondaryDepIn",&eProtonSecondaryDepIn.at(stop),"eProtonSecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("eProtonPrimaryDepOut",&eProtonPrimaryDepOut.at(stop),"eProtonPrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("eProtonSecondaryDepOut",&eProtonSecondaryDepOut.at(stop),"eProtonSecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("eNeutronPrimaryDepIn",&eNeutronPrimaryDepIn.at(stop),"eNeutronPrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("eNeutronSecondaryDepIn",&eNeutronSecondaryDepIn.at(stop),"eNeutronSecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("eNeutronPrimaryDepOut",&eNeutronPrimaryDepOut.at(stop),"eNeutronPrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("eNeutronSecondaryDepOut",&eNeutronSecondaryDepOut.at(stop),"eNeutronSecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("ePiCPrimaryDepIn",&ePiCPrimaryDepIn.at(stop),"ePiCPrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("ePiCSecondaryDepIn",&ePiCSecondaryDepIn.at(stop),"ePiCSecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("ePiCPrimaryDepOut",&ePiCPrimaryDepOut.at(stop),"ePiCPrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("ePiCSecondaryDepOut",&ePiCSecondaryDepOut.at(stop),"ePiCSecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("eGammaPrimaryDepIn",&eGammaPrimaryDepIn.at(stop),"eGammaPrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("eGammaSecondaryDepIn",&eGammaSecondaryDepIn.at(stop),"eGammaSecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("eGammaPrimaryDepOut",&eGammaPrimaryDepOut.at(stop),"eGammaPrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("eGammaSecondaryDepOut",&eGammaSecondaryDepOut.at(stop),"eGammaSecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("ePi0PrimaryDepIn",&ePi0PrimaryDepIn.at(stop),"ePi0PrimaryDepIn/D");
  stopsTrees.at(stop)->Branch("ePi0SecondaryDepIn",&ePi0SecondaryDepIn.at(stop),"ePi0SecondaryDepIn/D");
  stopsTrees.at(stop)->Branch("ePi0PrimaryDepOut",&ePi0PrimaryDepOut.at(stop),"ePi0PrimaryDepOut/D");
  stopsTrees.at(stop)->Branch("ePi0SecondaryDepOut",&ePi0SecondaryDepOut.at(stop),"ePi0SecondaryDepOut/D");

  stopsTrees.at(stop)->Branch("eReco",&eReco.at(stop),"eReco/D");

  stopsTrees.at(stop)->Branch("eHadTrueCharged",&eHadTrueCharged.at(stop),"eHadTrueCharged/D");
  stopsTrees.at(stop)->Branch("eHadTrueTotal",&eHadTrueTotal.at(stop),"eHadTrueTotal/D");
  stopsTrees.at(stop)->Branch("eMuTrue",&eMuTrue.at(stop),"eMuTrue/D");
  stopsTrees.at(stop)->Branch("eGammaTrue",&eGammaTrue.at(stop),"eGammaTrue/D");
  stopsTrees.at(stop)->Branch("pMuTrueX",&pMuTrueX.at(stop),"pMuTrueX/D");
  stopsTrees.at(stop)->Branch("pMuTrueY",&pMuTrueY.at(stop),"pMuTrueY/D");
  stopsTrees.at(stop)->Branch("pMuTrueZ",&pMuTrueZ.at(stop),"pMuTrueZ/D");

  stopsTrees.at(stop)->Branch("Q2True",&Q2True.at(stop),"Q2True/D"); 
  stopsTrees.at(stop)->Branch("yTrue",&yTrue.at(stop),"yTrue/D"); 
  stopsTrees.at(stop)->Branch("W_rest",&W_rest.at(stop),"W_rest/D"); 

  stopsTrees.at(stop)->Branch("nMu",&nMu.at(stop),"nMu/I");
  stopsTrees.at(stop)->Branch("nGamma",&nGamma.at(stop),"nGamma/I");
  stopsTrees.at(stop)->Branch("nPi0",&nPi0.at(stop),"nPi0/I");
  stopsTrees.at(stop)->Branch("nPiC",&nPiC.at(stop),"nPiC/I");
  stopsTrees.at(stop)->Branch("nProton",&nProton.at(stop),"nProton/I");
  stopsTrees.at(stop)->Branch("nNeutron",&nNeutron.at(stop),"nNeutron/I");

}

void DunePrismAnalyzer::InitVarsStop(){
  eventNum.push_back(0);
  Enu.push_back(0.);
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
  
  eMuPrimaryDep.push_back(0.);
  eMuSecondaryDep.push_back(0.);

  eHadPrimaryDepIn.push_back(0.);
  eHadSecondaryDepIn.push_back(0.);
  eHadPrimaryDepOut.push_back(0.);
  eHadSecondaryDepOut.push_back(0.);

  eProtonPrimaryDepIn.push_back(0.);
  eProtonSecondaryDepIn.push_back(0.);
  eProtonPrimaryDepOut.push_back(0.);
  eProtonSecondaryDepOut.push_back(0.);

  eNeutronPrimaryDepIn.push_back(0.);
  eNeutronSecondaryDepIn.push_back(0.);
  eNeutronPrimaryDepOut.push_back(0.);
  eNeutronSecondaryDepOut.push_back(0.);

  ePiCPrimaryDepIn.push_back(0.);
  ePiCSecondaryDepIn.push_back(0.);
  ePiCPrimaryDepOut.push_back(0.);
  ePiCSecondaryDepOut.push_back(0.);

  eGammaPrimaryDepIn.push_back(0.);
  eGammaSecondaryDepIn.push_back(0.);
  eGammaPrimaryDepOut.push_back(0.);
  eGammaSecondaryDepOut.push_back(0.);

  ePi0PrimaryDepIn.push_back(0.);
  ePi0SecondaryDepIn.push_back(0.);
  ePi0PrimaryDepOut.push_back(0.);
  ePi0SecondaryDepOut.push_back(0.);

  eReco.push_back(0.);

  eHadTrueCharged.push_back(0.);
  eHadTrueTotal.push_back(0.);
  eMuTrue.push_back(0.);
  eGammaTrue.push_back(0.);
  pMuTrueX.push_back(0.);
  pMuTrueY.push_back(0.);
  pMuTrueZ.push_back(0.);

  Q2True.push_back(0.);
  yTrue.push_back(0.);
  W_rest.push_back(0.);

  nMu.push_back(0.);
  nGamma.push_back(0.);
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
