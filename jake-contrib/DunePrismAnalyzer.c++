#include "DunePrismAnalyzer.h"

#ifndef DEBUG
 #define DEBUG
#endif



int main(int argc, char* argv[]){
  parse_args_xml(argc, argv);

//  gROOT->ProcessLine("#include <vector>");
//  gROOT->ProcessLine("#include <array>");

//  DunePrismAnalyzer dpa = DunePrismAnalyzer(inFileName, outFileName, detector, fiducialGap, offset);
  DunePrismAnalyzer dpa = DunePrismAnalyzer(inFileName, outFileName, detStops,nEntries);
  std::cout << "Checking dpa size: " << dpa.stopsTrees.size()<< std::endl;
  dpa.AnalyzeStops();
  dpa.FinalizeStops();
  return 0;
}


DunePrismAnalyzer::DunePrismAnalyzer(std::string inFileName, std::string outFileName, std::vector<DetectorStop> detStops, int N){
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

  //nEntries = inTree->GetEntries();
  nEntries = N;
  if (nEntries > inTree->GetEntries() || nEntries == -1) nEntries = inTree->GetEntries();
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

    InitDetector(i);

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
//    if(nuPID != 14) continue;

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
    nuPDG.at(stop_num) = nuPID;    
    vtx_X.at(stop_num) = EvtVtx[0];
    vtx_Y.at(stop_num) = EvtVtx[1];
    vtx_Z.at(stop_num) = EvtVtx[2];
    eventNum.at(stop_num) = ev; 

    double leptonBoundX[2] = {
      wallX[0] - FV.at(stop_num)[0],
      wallX[1] + FV.at(stop_num)[0]
    };
    double leptonBoundY = wallY + FV.at(stop_num)[1];
    double leptonBoundZ = wallZ + FV.at(stop_num)[1];
   

    eHadTrueCharged.at(stop_num) = 0.;
    eHadTrueTotal.at(stop_num) = 0.;
    eLepTrue.at(stop_num) = 0.;
    eGammaTrue.at(stop_num) = 0.;

    eProtonTrue.at(stop_num) = 0.;
    eNeutronTrue.at(stop_num) = 0.;
    ePi0True.at(stop_num) = 0.;
    ePiCTrue.at(stop_num) = 0.;


    pLepTrueX.at(stop_num) = 0.;
    pLepTrueY.at(stop_num) = 0.;
    pLepTrueZ.at(stop_num) = 0.;

    nLep.at(stop_num) = 0;
    nPi0.at(stop_num) = 0;
    nPiC.at(stop_num) = 0;
    nProton.at(stop_num) = 0;
    nNeutron.at(stop_num) = 0;
    nGamma.at(stop_num) = 0;

    std::map<int,int> chain;

    int nBindino = 0;
    int nOther = 0;

    FSHadrons.clear();

//    std::cout << "start of part loop" << std::endl;
    //Start of particle loop
    for(int ip = 0; ip < ni; ++ip){
      //std::cout << ip << std::endl;
      if (PIDi[ip] == 2000000101){
        nBindino++;
        continue;//Skip bindino
      }
      chain[ip + 1] = 0;       
      if((abs(PIDi[ip]) == 13 || abs(PIDi[ip]) == 11 || 
         PIDi[ip] == nuPID) && nLep.at(stop_num) < 1){
 //       std::cout << "Found " << PIDi[ip] << " " << ip + 1 << std::endl; 
        lepPDG.at(stop_num) = PIDi[ip];
        FSLepton = new DepoLepton(
          PIDi[ip], (ip + 1), ekini[ip] + mi[ip],
          leptonBoundX, leptonBoundY, leptonBoundZ); 

        eLepTrue.at(stop_num) = ekini[ip] + mi[ip];
        pLepTrueX.at(stop_num) = pxi[ip];
        pLepTrueY.at(stop_num) = pyi[ip];
        pLepTrueZ.at(stop_num) = pzi[ip];

        Q2True.at(stop_num) = fabs( pow( (ekina - eLepTrue.at(stop_num)), 2 ) 
                                 - pow( (pxi[ip] - pxa), 2 )
                                 - pow( (pyi[ip] - pya), 2 )
                                 - pow( (pzi[ip] - pza), 2 ) );

        yTrue.at(stop_num) = 1  - eLepTrue.at(stop_num)/ekina;
        double mN = .93827208;
        W_rest.at(stop_num) = sqrt(-Q2True.at(stop_num) + 2 * mN * (ekina - eLepTrue.at(stop_num)) + mN * mN);

        nLep.at(stop_num)++;

        if( abs(PIDi[ip]) == abs(nuPID) - 1 ) flagCC.at(stop_num) = 1;
        else if(PIDi[ip] == nuPID) flagCC.at(stop_num) = 0;

      }
      else if(PIDi[ip] == 111){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        ePi0True.at(stop_num) += ekini[ip] + mi[ip];
        nPi0.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip] + mi[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(abs(PIDi[ip]) == 211 ){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip] + mi[ip];
        ePiCTrue.at(stop_num) += ekini[ip] + mi[ip];
        nPiC.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip] + mi[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(PIDi[ip] == 2212){
//        std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueCharged.at(stop_num) += ekini[ip];
        eProtonTrue.at(stop_num) += ekini[ip];
        nProton.at(stop_num)++;
        DepoHadron * hadron = new DepoHadron(
          PIDi[ip], (ip + 1), ekini[ip],
          wallX, wallY, wallZ);
        FSHadrons[ip + 1] = hadron;
      }
      else if(PIDi[ip] == 2112){
 //       std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
        eHadTrueTotal.at(stop_num) += ekini[ip];
        eNeutronTrue.at(stop_num) += ekini[ip];
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
//    std::cout << "end of part loop" << std::endl;
    int nPrimary = ni - nBindino - nOther;
//    std::cout << "NHadrons: " << FSHadrons.size() << std::endl;
    if (nLep.at(stop_num) == 0) continue;

    eHadTrueTotal.at(stop_num) += eHadTrueCharged.at(stop_num); 
    ///////////////////////////////////

    flagExitBack.at(stop_num) = 0;
    flagExitFront.at(stop_num) = 0;
    flagExitY.at(stop_num) = 0;
    flagExitXHigh.at(stop_num) = 0;
    flagExitXLow.at(stop_num) = 0;
    flagLepContained.at(stop_num) = 0;
    flagNoEHadOut.at(stop_num) = 0;

    lepExitingPX.at(stop_num) = 0.;
    lepExitingPY.at(stop_num) = 0.;
    lepExitingPZ.at(stop_num) = 0.;

    eLepPrimaryDep.at(stop_num) = 0.;
    eLepSecondaryDep.at(stop_num) = 0.;

    eHadPrimaryDepIn.at(stop_num) = 0.;
    eHadSecondaryDepIn.at(stop_num) = 0.;
    eHadPrimaryDepOut.at(stop_num) = 0.;
    eHadSecondaryDepOut.at(stop_num) = 0.;

    for(int it = 0; it < xBins.at(stop_num).size() -1; ++it){
      for(int jt = 0; jt < yBins.at(stop_num).size() - 1; ++jt){
        for(int kt = 0; kt < zBins.at(stop_num).size() - 1; ++kt){
          eHadPrimaryDep[stop_num][it][jt][kt] = 0.;
          eProtonPrimaryDep[stop_num][it][jt][kt] = 0.;
          eNeutronPrimaryDep[stop_num][it][jt][kt] = 0.;
          ePiCPrimaryDep[stop_num][it][jt][kt] = 0.;
          ePi0PrimaryDep[stop_num][it][jt][kt] = 0.;
          eOtherPrimaryDep[stop_num][it][jt][kt] = 0.;

          eHadSecondaryDep[stop_num][it][jt][kt] = 0.;
          eProtonSecondaryDep[stop_num][it][jt][kt] = 0.;
          eNeutronSecondaryDep[stop_num][it][jt][kt] = 0.;
          ePiCSecondaryDep[stop_num][it][jt][kt] = 0.;
          ePi0SecondaryDep[stop_num][it][jt][kt] = 0.;
          eOtherSecondaryDep[stop_num][it][jt][kt] = 0.;
        }
      }
    }
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

    eOtherDepIn.at(stop_num) = 0.;
    eOtherDepOut.at(stop_num) = 0.;
    eTotalDep.at(stop_num) = 0.;
   
//    std::cout << "start steps" << std::endl;
    for(int i = 0; i < nstep; ++i){
      eTotalDep.at(stop_num) += edep[i];
      if(chain.find(track[i]) == chain.end()){//Not in chain
      
        if(chain.find(parid[i]) == chain.end()){//Error
          std::cout << "ERROR: PARENT NOT FOUND" << std::endl;
          break;
        }
        chain[track[i]] = parid[i];
      }
    
      if(chain[track[i]] == 0){//Primary
//        std::cout << "primary" << std::endl;
        if(PID[i] == lepPDG.at(stop_num)){//lepton
          if(FSLepton->trackID != track[i]){
            std::cout << "ERROR: Wrong lepton track" << std::endl;  
          }
          if(!FSLepton->flagLepContained)continue;
          FSLepton->eDepPrimary += edep[i];
          FSLepton->xf = xe[i];
          FSLepton->yf = ye[i];
          FSLepton->zf = ze[i];

          FSLepton->pxf = pxe[i];
          FSLepton->pyf = pye[i];
          FSLepton->pzf = pze[i];
          FSLepton->CheckContained();
        }
        else{//Hadron
          if(PID[i] == 2212 || PID[i] == 2112 || abs(PID[i]) == 211 ||
               PID[i] == 111  || PID[i] ==   22){

            if((xe[i] <= FSHadrons[track[i]]->xBound[0] && xe[i] >= FSHadrons[track[i]]->xBound[0] - FV.at(stop_num)[0]) || 
               (xe[i] >= FSHadrons[track[i]]->xBound[1] && xe[i] <= FSHadrons[track[i]]->xBound[1] + FV.at(stop_num)[0]) ||
               ( fabs(ye[i]) >= FSHadrons[track[i]]->yBound && fabs(ye[i]) <= FSHadrons[track[i]]->yBound + FV.at(stop_num)[1]) ||
               ( fabs(ze[i]) >= FSHadrons[track[i]]->zBound && fabs(ze[i]) <= FSHadrons[track[i]]->zBound + FV.at(stop_num)[2]) ){
              FSHadrons[track[i]]->eDepPrimaryOut += edep[i];
            }
            else{
//              std::cout << "in"<<std::endl;  
              FSHadrons[track[i]]->eDepPrimaryIn += edep[i];
            }         

            if( xe[i] >= xBins.at(stop_num).at(0) && xe[i] <= xBins.at(stop_num).at(xBins.at(stop_num).size() - 1) 
            &&  ye[i] >= yBins.at(stop_num).at(0) && ye[i] <= yBins.at(stop_num).at(3)
            &&  ze[i] >= zBins.at(stop_num).at(0) && ze[i] <= zBins.at(stop_num).at(3)){

              int bins[3] = {
                GetBinX(stop_num, xe[i]),
                GetBinY(stop_num, ye[i]),
                GetBinZ(stop_num, ze[i])};

              eHadPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
              switch (abs(PID[i])){
                case 2212:
                  eProtonPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case 2112:
                  eNeutronPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case  211:
                  ePiCPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case  111:
                  ePi0PrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                default:
                  eOtherPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
              }
            }

          }
          else{
            if((xe[i] <= wallX[0] && xe[i] >= wallX[0] - FV.at(stop_num)[0]) ||
               (xe[i] >= wallX[1] && xe[i] <= wallX[1] + FV.at(stop_num)[0]) ||
               ( fabs(ye[i]) >= wallY && fabs(ye[i]) <= wallY + FV.at(stop_num)[1]) ||
               ( fabs(ze[i]) >= wallZ && fabs(ze[i]) <= wallZ + FV.at(stop_num)[2]) ){
              eOtherDepOut.at(stop_num) += edep[i];
            }
            else{
              eOtherDepIn.at(stop_num) += edep[i];
            }
            if( xe[i] >= xBins.at(stop_num).at(0) && xe[i] <= xBins.at(stop_num).at(xBins.at(stop_num).size() - 1) 
            &&  ye[i] >= yBins.at(stop_num).at(0) && ye[i] <= yBins.at(stop_num).at(3)
            &&  ze[i] >= zBins.at(stop_num).at(0) && ze[i] <= zBins.at(stop_num).at(3)){

              int bins[3] = {
                GetBinX(stop_num, xe[i]),
                GetBinY(stop_num, ye[i]),
                GetBinZ(stop_num, ze[i])};

              eOtherPrimaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];              
            }
          }
        }
      }
      else{//Secondary
//        std::cout << "secondary" << std::endl;
        //Find ultimate parent
        int itChain = track[i];
        while (chain[itChain] != 0){
          itChain = chain[itChain];                
        }
        if(itChain == 1){//Lepton
          FSLepton->eDepSecondary += edep[i];
        }
        else{//Hadron

          if(PID[i] == 2212 || PID[i] == 2112 || abs(PID[i]) == 211 || PID[i] == 111  || PID[i] == 22){

            if(FSHadrons.find(itChain) == FSHadrons.end()){
//              std::cout << "ERROR: ULTIMATE HADRON NOT IN FS " << PIDi[itChain - 1] << " " << PID[i] << " " << edep[i] << std::endl;
              if((xe[i] <= wallX[0] && xe[i] >= wallX[0] - FV.at(stop_num)[0]) ||
                (xe[i] >= wallX[1] && xe[i] <= wallX[1] + FV.at(stop_num)[0]) ||
                ( fabs(ye[i]) >= wallY && fabs(ye[i]) <= wallY + FV.at(stop_num)[1]) ||
                ( fabs(ze[i]) >= wallZ && fabs(ze[i]) <= wallZ + FV.at(stop_num)[2]) ){
                eOtherDepOut.at(stop_num) += edep[i];
              }
              else{
                eOtherDepIn.at(stop_num) += edep[i];
              }
              continue;
            }
//            std::cout << "FOUND ULTIMATE " << FSHadrons[itChain]->PDG <<  " " << edep[i] << " " << PID[i] << std::endl; 
            if( (xe[i] <= FSHadrons[itChain]->xBound[0] && xe[i] >= FSHadrons[itChain]->xBound[0] - FV.at(stop_num)[0]) ||
                (xe[i] >= FSHadrons[itChain]->xBound[1] && xe[i] <= FSHadrons[itChain]->xBound[1] + FV.at(stop_num)[0]) ||
                (fabs(ye[i]) >= FSHadrons[itChain]->yBound && fabs(ye[i]) <= FSHadrons[itChain]->yBound + FV.at(stop_num)[1]) ||
                (fabs(ze[i]) >= FSHadrons[itChain]->zBound && fabs(ye[i]) <= FSHadrons[itChain]->zBound + FV.at(stop_num)[2]) ){
              
              FSHadrons[itChain]->eDepSecondaryOut += edep[i];
            }
            else{
              FSHadrons[itChain]->eDepSecondaryIn += edep[i];
            }
            if( xe[i] >= xBins.at(stop_num).at(0) && xe[i] <= xBins.at(stop_num).at(xBins.at(stop_num).size() - 1) 
            &&  ye[i] >= yBins.at(stop_num).at(0) && ye[i] <= yBins.at(stop_num).at(3)
            &&  ze[i] >= zBins.at(stop_num).at(0) && ze[i] <= zBins.at(stop_num).at(3)){

              int bins[3] = {
                GetBinX(stop_num, xe[i]),
                GetBinY(stop_num, ye[i]),
                GetBinZ(stop_num, ze[i])};

              eHadSecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
              switch (abs(FSHadrons[itChain]->PDG)){
                case 2212:
                  eProtonSecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case 2112:
                  eNeutronSecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case  211:
                  ePiCSecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;
                case  111:
                  ePi0SecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];
                  break;          
              }
            }
          }
          else{
            if( (xe[i] <= wallX[0] && xe[i] >= wallX[0] - FV.at(stop_num)[0]) ||
                (xe[i] >= wallX[1] && xe[i] <= wallX[1] + FV.at(stop_num)[0]) ||
                (fabs(ye[i]) >= wallY && fabs(ye[i]) <= wallY + FV.at(stop_num)[1]) ||
                (fabs(ze[i]) >= wallZ && fabs(ze[i]) <= wallZ + FV.at(stop_num)[2]) ){
              eOtherDepOut.at(stop_num) += edep[i];
            }
            else{
              eOtherDepIn.at(stop_num) += edep[i];
            }
            if( xe[i] >= xBins.at(stop_num).at(0) && xe[i] <= xBins.at(stop_num).at(xBins.at(stop_num).size() - 1) 
            &&  ye[i] >= yBins.at(stop_num).at(0) && ye[i] <= yBins.at(stop_num).at(3)
            &&  ze[i] >= zBins.at(stop_num).at(0) && ze[i] <= zBins.at(stop_num).at(3)){

              int bins[3] = {
                GetBinX(stop_num, xe[i]),
                GetBinY(stop_num, ye[i]),
                GetBinZ(stop_num, ze[i])};

              eOtherSecondaryDep[stop_num][bins[0]][bins[1]][bins[2]] += edep[i];              
            }
          }
        }
      }
    }
//    std::cout << "Ending steps " << std::endl;
    eLepPrimaryDep.at(stop_num) += FSLepton->eDepPrimary;
    eLepSecondaryDep.at(stop_num) += FSLepton->eDepSecondary;

    flagExitXLow.at(stop_num) = FSLepton->flagExitXLow;
    flagExitXHigh.at(stop_num) = FSLepton->flagExitXHigh;
    flagExitY.at(stop_num) = FSLepton->flagExitY;
    flagExitBack.at(stop_num) = FSLepton->flagExitBack;
    flagExitFront.at(stop_num) = FSLepton->flagExitFront;
    flagLepContained.at(stop_num) = FSLepton->flagLepContained;

    lepExitingPX.at(stop_num) = FSLepton->pxf;
    lepExitingPY.at(stop_num) = FSLepton->pyf;
    lepExitingPZ.at(stop_num) = FSLepton->pzf;   

    flagNoEHadOut.at(stop_num) = 0;
    
    std::map<int,DepoHadron*>::iterator itHad;
    int count = 0;
    for(itHad = FSHadrons.begin(); itHad != FSHadrons.end(); ++itHad){
      count++;

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

    if(count != FSHadrons.size()){std::cout << "Count: " << count << std::endl;    std::cout << "Nhads: " << FSHadrons.size() << std::endl;}

    eReco.at(stop_num) = eLepPrimaryDep.at(stop_num) + eLepSecondaryDep.at(stop_num) + eHadPrimaryDepIn.at(stop_num) + eHadSecondaryDepIn.at(stop_num) + eOtherDepIn.at(stop_num);


/*    for(int it = 0; it < 60; ++it){
      for(int jt = 0; jt < 3; ++jt){
        for(int kt = 0; kt < 3; ++kt){
          std::cout << eHadPrimaryDep[stop_num][it][jt][kt] << " ";
        }
      }
    }*/

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
  stopsTrees.at(stop)->Branch("nuPDG",&nuPDG.at(stop));
  stopsTrees.at(stop)->Branch("lepPDG",&lepPDG.at(stop));
  stopsTrees.at(stop)->Branch("flagCC",&flagCC.at(stop));
  stopsTrees.at(stop)->Branch("vtx_X",&vtx_X.at(stop),"vtx_X/D");
  stopsTrees.at(stop)->Branch("vtx_Y",&vtx_Y.at(stop),"vtx_Y/D");
  stopsTrees.at(stop)->Branch("vtx_Z",&vtx_Z.at(stop),"vtx_Z/D");

  stopsTrees.at(stop)->Branch("flagExitBack",&flagExitBack.at(stop));
  stopsTrees.at(stop)->Branch("flagExitFront",&flagExitFront.at(stop));
  stopsTrees.at(stop)->Branch("flagExitY",&flagExitY.at(stop));
  stopsTrees.at(stop)->Branch("flagExitXLow",&flagExitXLow.at(stop));
  stopsTrees.at(stop)->Branch("flagExitXHigh",&flagExitXHigh.at(stop));
  stopsTrees.at(stop)->Branch("flagNoEHadOut",&flagNoEHadOut.at(stop));
  stopsTrees.at(stop)->Branch("flagLepContained",&flagLepContained.at(stop));

  stopsTrees.at(stop)->Branch("lepExitingPX",&lepExitingPX.at(stop),"lepExitingPX/D");
  stopsTrees.at(stop)->Branch("lepExitingPY",&lepExitingPY.at(stop),"lepExitingPY/D");
  stopsTrees.at(stop)->Branch("lepExitingPZ",&lepExitingPZ.at(stop),"lepExitingPZ/D");

  stopsTrees.at(stop)->Branch("eLepPrimaryDep",&eLepPrimaryDep.at(stop),"eLepPrimaryDep/D");
  stopsTrees.at(stop)->Branch("eLepSecondaryDep",&eLepSecondaryDep.at(stop),"eLepSecondaryDep/D");

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

  stopsTrees.at(stop)->Branch("eOtherDepIn",&eOtherDepIn.at(stop),"eOtherDepIn/D");
  stopsTrees.at(stop)->Branch("eOtherDepOut",&eOtherDepOut.at(stop),"eOtherDepOut/D");
  stopsTrees.at(stop)->Branch("eTotalDep",&eTotalDep.at(stop),"eTotalDep/D");


  stopsTrees.at(stop)->Branch("eReco",&eReco.at(stop),"eReco/D");

  stopsTrees.at(stop)->Branch("eHadTrueCharged",&eHadTrueCharged.at(stop),"eHadTrueCharged/D");
  stopsTrees.at(stop)->Branch("eHadTrueTotal",&eHadTrueTotal.at(stop),"eHadTrueTotal/D");
  stopsTrees.at(stop)->Branch("eLepTrue",&eLepTrue.at(stop),"eLepTrue/D");

  stopsTrees.at(stop)->Branch("eProtonTrue",&eProtonTrue.at(stop),"eProtonTrue/D");
  stopsTrees.at(stop)->Branch("eNeutronTrue",&eNeutronTrue.at(stop),"eNeutronTrue/D");
  stopsTrees.at(stop)->Branch("ePiCTrue",&ePiCTrue.at(stop),"ePiCTrue/D");
  stopsTrees.at(stop)->Branch("ePi0True",&ePi0True.at(stop),"ePi0True/D");

  stopsTrees.at(stop)->Branch("eGammaTrue",&eGammaTrue.at(stop),"eGammaTrue/D");
  stopsTrees.at(stop)->Branch("pLepTrueX",&pLepTrueX.at(stop),"pLepTrueX/D");
  stopsTrees.at(stop)->Branch("pLepTrueY",&pLepTrueY.at(stop),"pLepTrueY/D");
  stopsTrees.at(stop)->Branch("pLepTrueZ",&pLepTrueZ.at(stop),"pLepTrueZ/D");

  stopsTrees.at(stop)->Branch("Q2True",&Q2True.at(stop),"Q2True/D"); 
  stopsTrees.at(stop)->Branch("yTrue",&yTrue.at(stop),"yTrue/D"); 
  stopsTrees.at(stop)->Branch("W_rest",&W_rest.at(stop),"W_rest/D"); 

  stopsTrees.at(stop)->Branch("nLep",&nLep.at(stop),"nLep/I");
  stopsTrees.at(stop)->Branch("nGamma",&nGamma.at(stop),"nGamma/I");
  stopsTrees.at(stop)->Branch("nPi0",&nPi0.at(stop),"nPi0/I");
  stopsTrees.at(stop)->Branch("nPiC",&nPiC.at(stop),"nPiC/I");
  stopsTrees.at(stop)->Branch("nProton",&nProton.at(stop),"nProton/I");
  stopsTrees.at(stop)->Branch("nNeutron",&nNeutron.at(stop),"nNeutron/I");

}

void DunePrismAnalyzer::InitDetector(int stop){
  std::cout << "Initializing detector " << stop << std::endl;
  double voxSize = 10.;//10cm

  double sizeX = detStops.at(stop).detectorSizeX*100.;
  double sizeY = detStops.at(stop).detectorSizeY*100.;
  double sizeZ = detStops.at(stop).detectorSizeZ*100.;

  double gapX = detStops.at(stop).fiducialGapX*100.;
  double gapY = detStops.at(stop).fiducialGapY*100.;
  double gapZ = detStops.at(stop).fiducialGapZ*100.;

  double shift = detStops.at(stop).shift*100.;    
  int nBinsX = floor(sizeX/voxSize);
 
  std::vector<double> tempBins;

//  double *tempX = new double[nBinsX];

  for(int ix = 0; ix < nBinsX + 1; ++ix){
    tempBins.push_back(shift - sizeX/2. + ix*voxSize);
    std::cout << tempBins.at(ix) << " ";
  }

  std::cout << std::endl;
  xBins.push_back(tempBins);
 
  std::vector<double> tempBinsY;
  tempBinsY.push_back(-1*sizeY/2.);
  tempBinsY.push_back(gapY - sizeY/2.);
  tempBinsY.push_back(sizeY/2. - gapY);
  tempBinsY.push_back(sizeY/2.);
  yBins.push_back(tempBinsY);

  std::vector<double> tempBinsZ;
  tempBinsZ.push_back(-1*sizeZ/2.);
  tempBinsZ.push_back(gapZ - sizeZ/2.);
  tempBinsZ.push_back(sizeZ/2. - gapZ);
  tempBinsZ.push_back(sizeZ/2.);
  zBins.push_back(tempBinsZ);

  
  for(int i = 0; i < xBins.at(stop).size() -1; ++i){
    for(int j = 0; j < yBins.at(stop).size() - 1; ++j){
      for(int k = 0; k < zBins.at(stop).size() - 1; ++k){
        eHadPrimaryDep[stop][i][j][k] = 0.;
        eProtonPrimaryDep[stop][i][j][k] = 0.;
        eNeutronPrimaryDep[stop][i][j][k] = 0.;
        ePiCPrimaryDep[stop][i][j][k] = 0.;
        ePi0PrimaryDep[stop][i][j][k] = 0.;       
        eOtherPrimaryDep[stop][i][j][k] = 0.;       

        eHadSecondaryDep[stop][i][j][k] = 0.;
        eProtonSecondaryDep[stop][i][j][k] = 0.;
        eNeutronSecondaryDep[stop][i][j][k] = 0.;
        ePiCSecondaryDep[stop][i][j][k] = 0.;
        ePi0SecondaryDep[stop][i][j][k] = 0.;
        eOtherSecondaryDep[stop][i][j][k] = 0.;
      }
    }
  }
  

  stopsTrees.at(stop)->Branch("eHadPrimaryDep",eHadPrimaryDep[stop],"eHadPrimaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eProtonPrimaryDep",eProtonPrimaryDep[stop],"eProtonPrimaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eNeutronPrimaryDep",eNeutronPrimaryDep[stop],"eNeutronPrimaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("ePiCPrimaryDep",ePiCPrimaryDep[stop],"ePiCPrimaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("ePi0PrimaryDep",ePi0PrimaryDep[stop],"ePi0PrimaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eOtherPrimaryDep",eOtherPrimaryDep[stop],"eOtherPrimaryDep[100][3][3]/D");

  stopsTrees.at(stop)->Branch("eHadSecondaryDep",eHadSecondaryDep[stop],"eHadSecondaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eProtonSecondaryDep",eProtonSecondaryDep[stop],"eProtonSecondaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eNeutronSecondaryDep",eNeutronSecondaryDep[stop],"eNeutronSecondaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("ePiCSecondaryDep",ePiCSecondaryDep[stop],"ePiCSecondaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("ePi0SecondaryDep",ePi0SecondaryDep[stop],"ePi0SecondaryDep[100][3][3]/D");
  stopsTrees.at(stop)->Branch("eOtherSecondaryDep",eOtherSecondaryDep[stop],"eOtherSecondaryDep[100][3][3]/D");
}

  
void DunePrismAnalyzer::InitVarsStop(){
  eventNum.push_back(0);
  Enu.push_back(0.);
  nuPDG.push_back(0);
  lepPDG.push_back(0);
  flagCC.push_back(0);
  vtx_X.push_back(0.);
  vtx_Y.push_back(0.);
  vtx_Z.push_back(0.);
  flagExitBack.push_back(0);
  flagExitFront.push_back(0);
  flagExitY.push_back(0);
  flagExitXHigh.push_back(0);
  flagExitXLow.push_back(0);
  flagNoEHadOut.push_back(0);
  flagLepContained.push_back(0);

  lepExitingPX.push_back(0.);
  lepExitingPY.push_back(0.);
  lepExitingPZ.push_back(0.);
  
  eLepPrimaryDep.push_back(0.);
  eLepSecondaryDep.push_back(0.);

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
  eOtherDepIn.push_back(0.);
  eOtherDepOut.push_back(0.);
  eTotalDep.push_back(0.);


  eHadTrueCharged.push_back(0.);
  eHadTrueTotal.push_back(0.);
  eLepTrue.push_back(0.);

  eProtonTrue.push_back(0.);
  eNeutronTrue.push_back(0.);
  ePiCTrue.push_back(0.);
  ePi0True.push_back(0.);

  eGammaTrue.push_back(0.);
  pLepTrueX.push_back(0.);
  pLepTrueY.push_back(0.);
  pLepTrueZ.push_back(0.);

  Q2True.push_back(0.);
  yTrue.push_back(0.);
  W_rest.push_back(0.);

  nLep.push_back(0.);
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

int DunePrismAnalyzer::GetBinX(int stop, double x){
  int xbin;
  for(int ix = 0; ix < xBins.at(stop).size() - 1; ++ix){
    xbin = ix;
    if(x <= xBins.at(stop).at(ix + 1)) break;
  }
  return xbin;
}

int DunePrismAnalyzer::GetBinY(int stop, double y){
  int ybin;
  if (y <= yBins.at(stop).at(1)){ybin = 0;}
  else if (y <= yBins.at(stop).at(2)){ybin=1;}
  else if (y <= yBins.at(stop).at(3)){ybin=2;}
  return ybin;
}

int DunePrismAnalyzer::GetBinZ(int stop, double z){
  int zbin;
  if (z <= zBins.at(stop).at(1)){zbin = 0;}
  else if (z <= zBins.at(stop).at(2)){zbin=1;}
  else if (z <= zBins.at(stop).at(3)){zbin=2;}
  return zbin;
}

int * DunePrismAnalyzer::GetBins(int stop, double x, double y, double z){

  int xbin, ybin, zbin;
  for(int ix = 0; ix < xBins.at(stop).size() - 1; ++ix){
    xbin = ix;
    if (x <= xBins.at(stop).at(ix + 1))break;
  }
                
  if (z <= zBins.at(stop).at(1)){zbin = 0;}
  else if (z <= zBins.at(stop).at(2)){zbin = 1;}
  else if (z <= zBins.at(stop).at(3)){zbin = 2;}

  if (y <= yBins.at(stop).at(1)){ybin = 0;}
  else if (y <= yBins.at(stop).at(2)){ybin=1;}
  else if (y <= yBins.at(stop).at(3)){ybin=2;}

  static int bins[3] = {xbin,ybin,zbin};
  return bins; 
}
