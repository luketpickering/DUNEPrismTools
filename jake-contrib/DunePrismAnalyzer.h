#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TXMLEngine.h"
#include "TH2.h"
#include "xml_parse.h"
#include "DepoParticle.h"
#include "TMap.h"
#include "TArray.h"
#include "TROOT.h"

#include <array>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>

//Dumb bullshit initialization from pyroot stuff 
const int maxInit = 100;
const int maxTracks = 100000;
const int maxNQ = 10000000;

class DunePrismAnalyzer{
  public:

    //Functions
    DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off,int N);
    DunePrismAnalyzer(std::string inFileName, std::string outFileName, std::vector<DetectorStop> detStops, DetectorStop fullDet, int N);
    virtual ~DunePrismAnalyzer(){}; 
    //void Analyze();
    void AnalyzeStops();
    void Finalize();
    void FinalizeStops();
    void SetBranches();
    void SetInBranches();
    void SetOutBranches( int stop);
    void InitVarsStop();
    void InitDetector();
    int * GetBins(double x, double y, double z);
    int GetBinX(double x);
    int GetBinY(double y);
    int GetBinZ(double z);
    ////////////////////

    TFile * fin;
    TFile * fout;

    TTree * inTree;
  
    DepoLepton * FSLepton;
    std::map<int,DepoHadron *> FSHadrons;
     
    std::vector<TTree*> stopsTrees;
    TTree * fullDetTree; 
    int nStops;

    std::vector< std::array<double,3> > dimension;
    std::vector< std::array<double,3> > FV;
    std::vector<double> shift;
  
    int nEntries;
    int ev;
    double ekina;
    double pxa;
    double pya;
    double pza;
    int nuPID;
    double EvtVtx[3];

    int nstep;
    int ni;
    int PIDi[maxInit];
    double pxi[maxInit];
    double pyi[maxInit];
    double pzi[maxInit];
    double ekini[maxInit];
    double mi[maxInit];

    int PID[maxNQ];
    int track[maxNQ];
    int parid[maxNQ];
    double xs[maxNQ];
    double xe[maxNQ];
    double ys[maxNQ];
    double ye[maxNQ];
    double zs[maxNQ];
    double ze[maxNQ];
    double pxs[maxNQ];
    double pxe[maxNQ];
    double pys[maxNQ];
    double pye[maxNQ];
    double pzs[maxNQ];
    double pze[maxNQ];
    double ekin[maxNQ];
    double edep[maxNQ];


    //Output
    std::vector<int> eventNum;
    std::vector<int> nuPDG;
    std::vector<int> lepPDG;
    std::vector<int> flagCC; 
    std::vector<double> Enu;
    std::vector<double> vtx_X;
    std::vector<double> vtx_Y;
    std::vector<double> vtx_Z;

    std::vector<int> flagExitBack;
    std::vector<int> flagExitFront;
    std::vector<int> flagExitY;
    std::vector<int> flagExitXHigh;
    std::vector<int> flagExitXLow;
    std::vector<int> flagNoEHadOut;
    std::vector<int> flagLepContained; 

    std::vector<double> lepExitingPX;
    std::vector<double> lepExitingPY;
    std::vector<double> lepExitingPZ;

    std::vector<double> eLepPrimaryDep;
    std::vector<double> eLepSecondaryDep;

    std::vector<double> eHadPrimaryDepIn;
    std::vector<double> eHadPrimaryDepOut;
    std::vector<double> eHadSecondaryDepIn;
    std::vector<double> eHadSecondaryDepOut;
    
    std::vector<double> eProtonPrimaryDepIn;
    std::vector<double> eProtonPrimaryDepOut;
    std::vector<double> eProtonSecondaryDepIn;
    std::vector<double> eProtonSecondaryDepOut;

    std::vector<double> eNeutronPrimaryDepIn;
    std::vector<double> eNeutronPrimaryDepOut;
    std::vector<double> eNeutronSecondaryDepIn;
    std::vector<double> eNeutronSecondaryDepOut;

    std::vector<double> ePiCPrimaryDepIn;
    std::vector<double> ePiCPrimaryDepOut;
    std::vector<double> ePiCSecondaryDepIn;
    std::vector<double> ePiCSecondaryDepOut;

    std::vector<double> eGammaPrimaryDepIn;
    std::vector<double> eGammaPrimaryDepOut;
    std::vector<double> eGammaSecondaryDepIn;
    std::vector<double> eGammaSecondaryDepOut;

    std::vector<double> ePi0PrimaryDepIn;
    std::vector<double> ePi0PrimaryDepOut;
    std::vector<double> ePi0SecondaryDepIn;
    std::vector<double> ePi0SecondaryDepOut;

    std::vector<double> eReco;
    std::vector<double> eOtherDepIn;
    std::vector<double> eOtherDepOut;
    std::vector<double> eTotalDep;

    //Positional Deposits
    double eLepPrimaryDepFull[4000][3][3];
    double eHadPrimaryDep[4000][3][3];//4000XSegmengs/9YZ
    double eProtonPrimaryDep[4000][3][3];
    double eNeutronPrimaryDep[4000][3][3];
    double ePiCPrimaryDep[4000][3][3];
    double ePi0PrimaryDep[4000][3][3];
    double eOtherPrimaryDep[4000][3][3];

    double eLepSecondaryDepFull[4000][3][3];
    double eHadSecondaryDep[4000][3][3];//4000XSegmengs/9YZ
    double eProtonSecondaryDep[4000][3][3];
    double eNeutronSecondaryDep[4000][3][3];
    double ePiCSecondaryDep[4000][3][3];
    double ePi0SecondaryDep[4000][3][3];
    double eOtherSecondaryDep[4000][3][3];
    std::vector<double> xBins;
    std::vector<double> yBins;
    std::vector<double> zBins;

    //FullDet metadata
    double EnuFull;
    int nuPDGFull;
    double vtx_X_full;
    double vtx_Y_full;
    double vtx_Z_full;   
    int eventNumFull;
    int lepPDGFull;
    int nLepFull;
    double eLepTrueFull;
    double pLepTrueXFull;
    double pLepTrueYFull;
    double pLepTrueZFull;
    double Q2TrueFull;
    double yTrueFull;
    double W_rest_full;
    int flagCCFull;
    double eHadTrueChargedFull;
    double eHadTrueTotalFull;
    double ePi0TrueFull;
    int nPi0Full;
    double ePiCTrueFull;
    int nPiCFull;
    double eProtonTrueFull;
    int nProtonFull;
    double eNeutronTrueFull;
    int nNeutronFull;
    double eGammaTrueFull;
    int nGammaFull;

    //////////////////////

    std::vector<double> eHadTrueCharged;
    std::vector<double> eHadTrueTotal;
    std::vector<double> eLepTrue;

    std::vector<double> eProtonTrue;
    std::vector<double> eNeutronTrue;
    std::vector<double> ePiCTrue;
    std::vector<double> ePi0True;

    std::vector<double> eGammaTrue;
    std::vector<double> pLepTrueX;
    std::vector<double> pLepTrueY;
    std::vector<double> pLepTrueZ;

    std::vector<double> Q2True;
    std::vector<double> yTrue;
    std::vector<double> W_rest;

    std::vector<int> nLep;
    std::vector<int> nGamma;
    std::vector<int> nPi0;
    std::vector<int> nPiC;
    std::vector<int> nProton;
    std::vector<int> nNeutron;

  private:
};




void parse_args(int, char ** );

double detector[3]; 
double fiducialGap[3];
double offset; 
std::string inFileName;
std::string outFileName;
int nEntries;
std::vector<DetectorStop> detStops;
DetectorStop fullDet;


void parse_args(int argc, char* argv[]){
  for (int i = 1; i < argc; ++i){

    std::string flag = argv[i];

    if(flag == "-i") inFileName = argv[i+1];

    else if(flag == "-o") outFileName = argv[i+1];
    else if(flag.find("--det") != std::string::npos){ 
      if(flag.back() == 'X') detector[0] = atof(argv[i+1]);
      if(flag.back() == 'Y') detector[1] = atof(argv[i+1]);
      if(flag.back() == 'Z') detector[2] = atof(argv[i+1]);
    }

    else if(flag.find("--fid") != std::string::npos){
      if(flag.back() == 'X') fiducialGap[0] = atof(argv[i+1]);
      if(flag.back() == 'Y') fiducialGap[1] = atof(argv[i+1]);
      if(flag.back() == 'Z') fiducialGap[2] = atof(argv[i+1]);
    }

    else if(flag == "--offset") offset = atof(argv[i+1]);
  }
}

void parse_args_xml(int argc, char* argv[]){
  for (int i = 1; i < argc; ++i){

    std::string flag = argv[i];

    if(flag == "-i") inFileName = argv[i+1];

    else if(flag == "-o") outFileName = argv[i+1];

    else if(flag == "-x"){
      std::cout << argv[i+1] << std::endl;
      detStops = ReadDetectorStopConfig(argv[i+1]);
      fullDet = GetFullDetectorConfig(argv[i+1]);
      std::cout << "size "<<  detStops.size() << std::endl;
    }
    else if(flag == "-n"){
      nEntries = atoi(argv[i+1]);
    }
    
  }
}



