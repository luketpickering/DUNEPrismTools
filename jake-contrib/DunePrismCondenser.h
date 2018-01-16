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

class DunePrismCondenser{
  public:

    //Functions
    DunePrismCondenser(std::string inFileName, std::string outFileName, DetectorStop * fullDet, int N);
    virtual ~DunePrismCondenser(){}; 
    //void Analyze();
    void Condense();
    void Finalize();
    void SetInBranches();
    void InitDetector();
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
    double te[maxNQ];
    double pxs[maxNQ];
    double pxe[maxNQ];
    double pys[maxNQ];
    double pye[maxNQ];
    double pzs[maxNQ];
    double pze[maxNQ];
    double ekin[maxNQ];
    double edep[maxNQ];



    //Output
    //Positional Deposits
    double eLepPrimaryDepFull[400][3][3];
    double eHadPrimaryDep[400][3][3];//400XSegmengs/9YZ
    double eProtonPrimaryDep[400][3][3];
    double eNeutronPrimaryDep[400][3][3];
    double ePiCPrimaryDep[400][3][3];
    double ePi0PrimaryDep[400][3][3];
    double eOtherPrimaryDep[400][3][3];

    double eLepSecondaryDepFull[400][3][3];
    double eHadSecondaryDep[400][3][3];//400XSegmengs/9YZ
    double eProtonSecondaryDep[400][3][3];
    double eNeutronSecondaryDep[400][3][3];
    double ePiCSecondaryDep[400][3][3];
    double ePi0SecondaryDep[400][3][3];
    double eOtherSecondaryDep[400][3][3];
    double neutronPrimaryTimeDep[400][3][3];
    double neutronSecondaryTimeDep[400][3][3];
    std::vector<double> xBins;
    std::vector<double> yBins;
    std::vector<double> zBins;

    //FullDet metadata
    double lepExitingPosX;
    double lepExitingPosY;
    double lepExitingPosZ;
    double lepExitingMomX;
    double lepExitingMomY;
    double lepExitingMomZ;
    bool flagLepExit;
    bool flagLepExitBack;
    bool flagLepExitFront;
    bool flagLepExitY;
    /*double lepTrackX[1000];
    double lepTrackY[1000];
    double lepTrackZ[1000];
    double lepTrackMomX[1000];
    double lepTrackMomY[1000];
    double lepTrackMomZ[1000];*/
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


  private:
};




void parse_args(int, char ** );

double detector[3]; 
double fiducialGap[3];
double offset; 
std::string inFileName;
std::string outFileName;
int nEntries;
DetectorStop * fullDet;


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
      fullDet = GetFullDetectorConfig(argv[i+1]);
      std::cout << "fullDet: " << fullDet << std::endl;
    }
    else if(flag == "-n"){
      nEntries = atoi(argv[i+1]);
    }
    
  }
}



