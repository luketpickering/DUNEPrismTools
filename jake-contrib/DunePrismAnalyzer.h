#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TXMLEngine.h"
#include "TH2.h"
#include "xml_parse.h"

#include <array>
#include <iostream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

//Dumb bullshit initialization from pyroot stuff 
const int maxInit = 100;
const int maxTracks = 100000;
const int maxNQ = 10000000;

class DunePrismAnalyzer{
  public:

    //Functions
    DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off);
    DunePrismAnalyzer(std::string inFileName, std::string outFileName, std::vector<DetectorStop> detStops);
    virtual ~DunePrismAnalyzer(){}; 
    //void Analyze();
    void AnalyzeStops();
    void Finalize();
    void FinalizeStops();
    void SetBranches();
    void SetInBranches();
    void SetOutBranches( int stop);
    void InitVarsStop();
    ////////////////////

    TFile * fin;
    TFile * fout;

    TTree * inTree;
    TTree * eventTree;
    TTree * gOutTree;
   
    std::vector<TTree*> stopsTrees;
    int nStops;

    std::vector<int> trackPi0;//trackIDs of pi0s
    std::map<int,int> trackGamma;//< trackID -> parid > of gammas
    std::vector<int> trackPiC;//trackIDs of pics
    std::vector<int> trackProton;//trackIDs of protons
    std::vector<int> trackMu;//trackIDs of muons
    std::vector< std::pair<int,double> > trackE; //(parid,edep)
    std::vector< std::pair<int,double> > trackEIn; //(parid,edep)
    std::vector< std::pair<int,double> > trackEOut; //(parid,edep)

/*    double dimension[3];
    double FV[3];
    double shift;*/

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

    std::vector<double> mu_x;
    std::vector<double> mu_y;
    std::vector<double> mu_z;
    std::vector<double> mu_px;
    std::vector<double> mu_py;
    std::vector<double> mu_pz;

    //Output
    std::vector<int> eventNum;
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
    std::vector<int> flagMuContained; 

    std::vector<double> muExitingPX;
    std::vector<double> muExitingPY;
    std::vector<double> muExitingPZ;

    std::vector<double> eHadOutDep;
    std::vector<double> eHadInDep;
    std::vector<double> eHadTotalDep;

    std::vector<double> ePi0TotalDep;
    std::vector<double> ePi0InDep;
    std::vector<double> ePi0OutDep;

    std::vector<double> ePiCTotalDep;
    std::vector<double> ePiCInDep;
    std::vector<double> ePiCOutDep;

    std::vector<double> eProtonTotalDep;
    std::vector<double> eProtonInDep;
    std::vector<double> eProtonOutDep;

    std::vector<double> ePiCEMTotalDep;
    std::vector<double> ePiCEMInDep;
    std::vector<double> ePiCEMOutDep;

    std::vector<double> eProtonEMTotalDep;
    std::vector<double> eProtonEMInDep;
    std::vector<double> eProtonEMOutDep;

    std::vector<double> eMuEMTotalDep;
    std::vector<double> eMuEMInDep;
    std::vector<double> eMuEMOutDep;

    std::vector<double> eResidualEMTotalDep;
    std::vector<double> eResidualEMInDep;
    std::vector<double> eResidualEMOutDep;

    std::vector<double> eMuDep;
    std::vector<double> eMuTotalDep;
    std::vector<double> eMuSecondaryDep;

    std::vector<double> eReco;

    std::vector<double> eHadTrueCharged;
    std::vector<double> eHadTrueTotal;
    std::vector<double> eMuTrue;
    std::vector<double> pMuTrueX;
    std::vector<double> pMuTrueY;
    std::vector<double> pMuTrueZ;

    std::vector<double> Q2True;
    std::vector<double> yTrue;
    std::vector<double> W_rest;

    std::vector<int> nMu;
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

std::vector<DetectorStop> detStops;


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
      std::cout << "size "<<  detStops.size() << std::endl;
    }
    
  }
}



