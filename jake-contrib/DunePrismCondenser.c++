#include "DunePrismCondenser.h"

#ifndef DEBUG
 #define DEBUG
#endif



int main(int argc, char* argv[]){
  parse_args_xml(argc, argv);

  DunePrismCondenser dpc = DunePrismCondenser(inFileName, outFileName,  fullDet,nEntries);
  dpc.Condense();
  dpc.Finalize();
  return 0;
}


DunePrismCondenser::DunePrismCondenser(std::string inFileName, std::string outFileName, DetectorStop * fullDet, int N){
 #ifdef DEBUG
  std::cout << inFileName << std::endl;
  std::cout << outFileName << std::endl;
  #endif

  fin = new TFile(inFileName.c_str(),"READ");
  inTree = (TTree*)fin->Get("argon"); 

  if(!inTree) std::cout << "Can't get argon tree" << std::endl;
  #ifdef DEBUG
  else std::cout<<"Got argon tree" << std::endl;
  #endif

  nEntries = N;
  if (nEntries > inTree->GetEntries() || nEntries == -1) nEntries = inTree->GetEntries();
  #ifdef DEBUG
  std::cout << nEntries << std::endl; 
  #endif

 
  fout = new TFile(outFileName.c_str(),"RECREATE");

  if(fullDet != NULL){
    fullDetTree = new TTree("fullDetTree","Positional deposits over full detector");
    InitDetector();
  }

  SetInBranches();
}

void DunePrismCondenser::Condense(){
  for (int ie = 0; ie < nEntries; ++ie){
    if(!(ie%100))std::cout<<"Entry " << ie<<std::endl;
    inTree->GetEntry(ie);
//    if(nuPID != 14) continue;


    double wallX[2] = {0.,0.};
    double wallY = 0.;
    double wallZ = 0.;
 

      //Start of full detector info 
      EnuFull = ekina;
      nuPDGFull = nuPID;    
      vtx_X_full = EvtVtx[0];
      vtx_Y_full = EvtVtx[1];
      vtx_Z_full = EvtVtx[2];
      eventNumFull = ev;

      if(vtx_X_full < (fullDet->shift - fullDet->detectorSizeX*.5)*100. ||
         vtx_X_full > (fullDet->shift + fullDet->detectorSizeX*.5)*100. ||
         fabs(vtx_Y_full) > fullDet->detectorSizeY*.5*100 ||
         fabs(vtx_Z_full) > fullDet->detectorSizeZ*.5*100){
           std::cout<<"Skipping "<<ie<<std::endl << "\t" <<
           vtx_X_full << " " << vtx_Y_full << " " << vtx_Z_full<<std::endl << "\t" << 
           (fullDet->shift - fullDet->detectorSizeX*.5)*100. << " " << (fullDet->shift + fullDet->detectorSizeX*.5)*100. << " " << fullDet->detectorSizeY*.5*100. << " " << fullDet->detectorSizeZ*.5*100. << std::endl;
           continue;
         }

        nLepFull = 0;
        lepPDGFull = 0;
        eLepTrueFull = 0.;
        pLepTrueXFull = 0.;
        pLepTrueYFull = 0.;
        pLepTrueZFull = 0.;
        Q2TrueFull = 0.;
        yTrueFull = 0.;
        W_rest_full = 0.;
        eHadTrueChargedFull = 0.;
        eHadTrueTotalFull = 0.;
        ePi0TrueFull = 0.;
        nPi0Full = 0;
        ePiCTrueFull = 0.;
        nPiCFull = 0;
        eProtonTrueFull = 0.;
        nProtonFull = 0;
        eNeutronTrueFull = 0.;
        nNeutronFull = 0;
        eGammaTrueFull = 0.;
        nGammaFull = 0;  

        std::map<int,int> chain;

        int nBindino = 0;
        int nOther = 0;

        FSHadrons.clear();
 
        //std::cout << "start of part loop" << std::endl;
        //Start of particle loop
        for(int ip = 0; ip < ni; ++ip){
          //std::cout << ip <<" " <<PIDi[ip]<< std::endl;
          if (PIDi[ip] == 2000000101){
            nBindino++;
            continue;//Skip bindino
          }
          chain[ip + 1] = 0;       
          if((abs(PIDi[ip]) == 13 || abs(PIDi[ip]) == 11 || 
             PIDi[ip] == nuPID) && nLepFull < 1){
            //std::cout << "Found " << PIDi[ip] << " " << ip + 1 << std::endl; 
            lepPDGFull = PIDi[ip];

            eLepTrueFull = ekini[ip] + mi[ip];
            pLepTrueXFull = pxi[ip];
            pLepTrueYFull = pyi[ip];
            pLepTrueZFull = pzi[ip];

            Q2TrueFull = fabs( pow( (ekina - eLepTrueFull), 2 ) 
                                     - pow( (pxi[ip] - pxa), 2 )
                                     - pow( (pyi[ip] - pya), 2 )
                                     - pow( (pzi[ip] - pza), 2 ) );

            yTrueFull = 1  - eLepTrueFull/ekina;
            double mN = .93827208;
            W_rest_full = sqrt(-Q2TrueFull + 2 * mN * (ekina - eLepTrueFull) + mN * mN);

            nLepFull++;

            if( abs(PIDi[ip]) == abs(nuPID) - 1 ) flagCCFull = 1;
            else if(PIDi[ip] == nuPID) flagCCFull = 0;

          }
          else if(PIDi[ip] == 111){
           // std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
            eHadTrueChargedFull += ekini[ip] + mi[ip];
            ePi0TrueFull += ekini[ip] + mi[ip];
            nPi0Full++;
            DepoHadron * hadron = new DepoHadron(
              PIDi[ip], (ip + 1), ekini[ip] + mi[ip],
              wallX, wallY, wallZ);
            FSHadrons[ip + 1] = hadron;
          }
          else if(abs(PIDi[ip]) == 211 ){
           // std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
            eHadTrueChargedFull += ekini[ip] + mi[ip];
            ePiCTrueFull += ekini[ip] + mi[ip];
            nPiCFull++;
            DepoHadron * hadron = new DepoHadron(
              PIDi[ip], (ip + 1), ekini[ip] + mi[ip],
              wallX, wallY, wallZ);
            FSHadrons[ip + 1] = hadron;
          }
          else if(PIDi[ip] == 2212){
           // std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
            eHadTrueChargedFull += ekini[ip];
            eProtonTrueFull += ekini[ip];
            nProtonFull++;
            DepoHadron * hadron = new DepoHadron(
              PIDi[ip], (ip + 1), ekini[ip],
              wallX, wallY, wallZ);
            FSHadrons[ip + 1] = hadron;
          }
          else if(PIDi[ip] == 2112){
         // std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
            eHadTrueTotalFull += ekini[ip];
            eNeutronTrueFull += ekini[ip];
            nNeutronFull++;
            DepoHadron * hadron = new DepoHadron(
              PIDi[ip], (ip + 1), ekini[ip],
              wallX, wallY, wallZ);
            FSHadrons[ip + 1] = hadron;
          }
          else if(PIDi[ip] == 22){
           // std::cout << "Found " << PIDi[ip] <<  " " << ip + 1 << std::endl; 
            eGammaTrueFull += ekini[ip];
            nGammaFull++;
            DepoHadron * hadron = new DepoHadron(
              PIDi[ip], (ip + 1), ekini[ip],
              wallX, wallY, wallZ);
            FSHadrons[ip + 1] = hadron;
            //Add it to hadrons because I'm lazy 
            //to make a new thing
          }
          else{
//            std::cout << "Found Other! " << PIDi[ip] <<  " " << ip + 1 << std::endl;
            nOther++;
          }
        }

     // std::cout << "Resetting bins\n";
      for(int it = 0; it < xBins.size() - 1; ++it){
        for(int jt = 0; jt < yBins.size() - 1; ++jt){
          for(int kt = 0; kt < zBins.size() - 1; ++kt){
            eLepPrimaryDepFull[it][jt][kt] = 0.;
            eHadPrimaryDep[it][jt][kt] = 0.;
            eProtonPrimaryDep[it][jt][kt] = 0.;
            eNeutronPrimaryDep[it][jt][kt] = 0.;
            ePiCPrimaryDep[it][jt][kt] = 0.;
            ePi0PrimaryDep[it][jt][kt] = 0.;
            eOtherPrimaryDep[it][jt][kt] = 0.;

            eLepSecondaryDepFull[it][jt][kt] = 0.;
            eHadSecondaryDep[it][jt][kt] = 0.;
            eProtonSecondaryDep[it][jt][kt] = 0.;
            eNeutronSecondaryDep[it][jt][kt] = 0.;
            ePiCSecondaryDep[it][jt][kt] = 0.;
            ePi0SecondaryDep[it][jt][kt] = 0.;
            eOtherSecondaryDep[it][jt][kt] = 0.;
          }
        }
      }

     // std::cout << "Steps\n";
      for(int i = 0; i < nstep; ++i){
        if(chain.find(track[i]) == chain.end()){//Not in chain
        
          if(chain.find(parid[i]) == chain.end()){//Error
            std::cout << "ERROR: PARENT NOT FOUND" << std::endl;
            break;
          }
          chain[track[i]] = parid[i];
        }
      
        if(chain[track[i]] == 0){//Primary
        //  std::cout << "primary" << std::endl;
          if(PID[i] == lepPDGFull){//lepton
        //    std::cout << "lepton" << std::endl;
            if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
               &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
               &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){  
               int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};
               eLepPrimaryDepFull[bins[0]][bins[1]][bins[2]] += edep[i];
               std::cout << ev << " " << bins[0] << " " << bins[1] << " " << bins[2] << " "<< edep[i]<<"\n";
            }
          }
          else{//Hadron
            if(PID[i] == 2212 || PID[i] == 2112 || abs(PID[i]) == 211 ||
                 PID[i] == 111  || PID[i] ==   22){
        //      std::cout << "Hadron" << std::endl;


              if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
              &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
              &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){

                int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};

                eHadPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                switch (abs(PID[i])){
                  case 2212:
                    eProtonPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case 2112:
                    eNeutronPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case  211:
                    ePiCPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case  111:
                    ePi0PrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  default:
                    eOtherPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                }
              }

            }
            else{
        //      std::cout << "Other" << std::endl;

              if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
              &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
              &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){

                int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};

                eOtherPrimaryDep[bins[0]][bins[1]][bins[2]] += edep[i];              
              }
            }
          }
        }
        else{//Secondary
         // std::cout << "secondary" << std::endl;
          //Find ultimate parent
          int itChain = track[i];
          while (chain[itChain] != 0){
            itChain = chain[itChain];                
          }
          if(itChain == 1){//Lepton
            //std::cout << "Lepton" << std::endl;
            if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
               &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
               &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){  
               int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};
               eLepSecondaryDepFull[bins[0]][bins[1]][bins[2]] += edep[i];
            }
          }
          else{//Hadron

            if( FSHadrons.find(itChain) != FSHadrons.end() ){ 
             // std::cout << "Hadron" << std::endl;

              if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
              &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
              &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){

                int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};

                eHadSecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
               // std::cout << itChain << " " ;
                //std::cout <<abs(FSHadrons[itChain]->PDG) << std::endl;
                switch (abs(FSHadrons[itChain]->PDG)){
                  case 2212:
                    eProtonSecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case 2112:
                    eNeutronSecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case  211:
                    ePiCSecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;
                  case  111:
                    ePi0SecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];
                    break;          
                }
              }
            }
            else{
          // std::cout << "Other" << std::endl;

              if( xe[i] >= xBins.at(0) && xe[i] <= xBins.at(xBins.size() - 1) 
              &&  ye[i] >= yBins.at(0) && ye[i] <= yBins.at(3)
              &&  ze[i] >= zBins.at(0) && ze[i] <= zBins.at(3)){

                int bins[3] = {
                  GetBinX(xe[i]),
                  GetBinY(ye[i]),
                  GetBinZ(ze[i])};

                eOtherSecondaryDep[bins[0]][bins[1]][bins[2]] += edep[i];              
              }

            }
          }
        }
      }
      //std::cout << "end Steps\n";
      fullDetTree->Fill();
    
  }
}


void DunePrismCondenser::SetInBranches(){

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



void DunePrismCondenser::InitDetector(){
  std::cout << "Initializing detector " << std::endl;
  double voxSize = 10.;//10cm

  double sizeX = fullDet->detectorSizeX*100.;
  double sizeY = fullDet->detectorSizeY*100.;
  double sizeZ = fullDet->detectorSizeZ*100.;

  double gapX = fullDet->fiducialGapX*100.;
  double gapY = fullDet->fiducialGapY*100.;
  double gapZ = fullDet->fiducialGapZ*100.;

  double shift = fullDet->shift*100.;    
  int nBinsX = floor(sizeX/voxSize);
 
  std::vector<double> tempBins;

//  double *tempX = new double[nBinsX];

  for(int ix = 0; ix < nBinsX + 1; ++ix){
    xBins.push_back(shift - sizeX/2. + ix*voxSize);
    std::cout << xBins.at(ix) << " ";
  }

  std::cout << std::endl;
 
  yBins.push_back(-1*sizeY/2.);
  yBins.push_back(gapY - sizeY/2.);
  yBins.push_back(sizeY/2. - gapY);
  yBins.push_back(sizeY/2.);

  zBins.push_back(-1*sizeZ/2.);
  zBins.push_back(gapZ - sizeZ/2.);
  zBins.push_back(sizeZ/2. - gapZ);
  zBins.push_back(sizeZ/2.);

  
  for(int i = 0; i < xBins.size() -1; ++i){
    for(int j = 0; j < yBins.size() - 1; ++j){
      for(int k = 0; k < zBins.size() - 1; ++k){
        eLepPrimaryDepFull[i][j][k] = 0.;
        eHadPrimaryDep[i][j][k] = 0.;
        eProtonPrimaryDep[i][j][k] = 0.;
        eNeutronPrimaryDep[i][j][k] = 0.;
        ePiCPrimaryDep[i][j][k] = 0.;
        ePi0PrimaryDep[i][j][k] = 0.;       
        eOtherPrimaryDep[i][j][k] = 0.;       

        eLepSecondaryDepFull[i][j][k] = 0.;
        eHadSecondaryDep[i][j][k] = 0.;
        eProtonSecondaryDep[i][j][k] = 0.;
        eNeutronSecondaryDep[i][j][k] = 0.;
        ePiCSecondaryDep[i][j][k] = 0.;
        ePi0SecondaryDep[i][j][k] = 0.;
        eOtherSecondaryDep[i][j][k] = 0.;
      }
    }
  }
  
  fullDetTree->Branch("Enu",&EnuFull);
  fullDetTree->Branch("nuPDG",&nuPDGFull);
  fullDetTree->Branch("eventNum",&eventNumFull);
  fullDetTree->Branch("vtx_X",&vtx_X_full);
  fullDetTree->Branch("vtx_Y",&vtx_Y_full);
  fullDetTree->Branch("vtx_Z",&vtx_Z_full);
  fullDetTree->Branch("lepPDG",&lepPDGFull);
  fullDetTree->Branch("nLep",&nLepFull);
  fullDetTree->Branch("eLepTrue",&eLepTrueFull);
  fullDetTree->Branch("pLepTrueX",&pLepTrueXFull);
  fullDetTree->Branch("pLepTrueY",&pLepTrueYFull);
  fullDetTree->Branch("pLepTrueZ",&pLepTrueZFull);
  fullDetTree->Branch("Q2True",&Q2TrueFull);
  fullDetTree->Branch("yTrue",&yTrueFull);
  fullDetTree->Branch("W_rest",&W_rest_full);
  fullDetTree->Branch("flagCC",&flagCCFull);
  fullDetTree->Branch("eHadTrueCharged",&eHadTrueChargedFull);
  fullDetTree->Branch("eHadTrueTotal",&eHadTrueTotalFull);
  fullDetTree->Branch("nPi0",&nPi0Full);
  fullDetTree->Branch("ePi0True",&ePi0TrueFull);
  fullDetTree->Branch("nPiC",&nPiCFull);
  fullDetTree->Branch("ePiCTrue",&ePiCTrueFull);
  fullDetTree->Branch("nProton",&nProtonFull);
  fullDetTree->Branch("eProtonTrue",&eProtonTrueFull);
  fullDetTree->Branch("nNeutron",&nNeutronFull);
  fullDetTree->Branch("eNeutronTrue",&eNeutronTrueFull);
  fullDetTree->Branch("nGamma",&nGammaFull);
  fullDetTree->Branch("eGammaTrue",&eGammaTrueFull);
 
  fullDetTree->Branch("eLepPrimaryDepFull",&eLepPrimaryDepFull,"eLepPrimaryDep[400][3][3]/D");
  fullDetTree->Branch("eHadPrimaryDep",&eHadPrimaryDep,"eHadPrimaryDep[400][3][3]/D");
  fullDetTree->Branch("eProtonPrimaryDep",&eProtonPrimaryDep,"eProtonPrimaryDep[400][3][3]/D");
  fullDetTree->Branch("eNeutronPrimaryDep",&eNeutronPrimaryDep,"eNeutronPrimaryDep[400][3][3]/D");
  fullDetTree->Branch("ePiCPrimaryDep",&ePiCPrimaryDep,"ePiCPrimaryDep[400][3][3]/D");
  fullDetTree->Branch("ePi0PrimaryDep",&ePi0PrimaryDep,"ePi0PrimaryDep[400][3][3]/D");
  fullDetTree->Branch("eOtherPrimaryDep",&eOtherPrimaryDep,"eOtherPrimaryDep[400][3][3]/D");

  fullDetTree->Branch("eLepSecondaryDepFull",&eLepSecondaryDepFull,"eLepSecondaryDep[400][3][3]/D");
  fullDetTree->Branch("eHadSecondaryDep",&eHadSecondaryDep,"eHadSecondaryDep[400][3][3]/D");
  fullDetTree->Branch("eProtonSecondaryDep",&eProtonSecondaryDep,"eProtonSecondaryDep[400][3][3]/D");
  fullDetTree->Branch("eNeutronSecondaryDep",&eNeutronSecondaryDep,"eNeutronSecondaryDep[400][3][3]/D");
  fullDetTree->Branch("ePiCSecondaryDep",&ePiCSecondaryDep,"ePiCSecondaryDep[400][3][3]/D");
  fullDetTree->Branch("ePi0SecondaryDep",&ePi0SecondaryDep,"ePi0SecondaryDep[400][3][3]/D");
  fullDetTree->Branch("eOtherSecondaryDep",&eOtherSecondaryDep,"eOtherSecondaryDep[400][3][3]/D");
}

  


void DunePrismCondenser::Finalize(){
  fout->cd();
  if(fullDet != NULL){
    fullDetTree->Write(0,TObject::kOverwrite);
  }
  fout->Close();
  fin->Close();

}

int DunePrismCondenser::GetBinX(double x){
  int xbin;
  for(int ix = 0; ix < xBins.size() - 1; ++ix){
    xbin = ix;
    if(x <= xBins.at(ix + 1)) break;
  }
  return xbin;
}

int DunePrismCondenser::GetBinY(double y){
  int ybin;
  if (y <= yBins.at(1)){ybin = 0;}
  else if (y <= yBins.at(2)){ybin=1;}
  else if (y <= yBins.at(3)){ybin=2;}
  return ybin;
}

int DunePrismCondenser::GetBinZ(double z){
  int zbin;
  if (z <= zBins.at(1)){zbin = 0;}
  else if (z <= zBins.at(2)){zbin=1;}
  else if (z <= zBins.at(3)){zbin=2;}
  return zbin;
}


