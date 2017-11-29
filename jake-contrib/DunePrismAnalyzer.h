#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <array>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

//Dumb bullshit initialization from pyroot stuff 
const int maxInit = 100;
const int maxTracks = 100000;
const int maxNQ = 10000000;

class DunePrismAnalyzer{
  public:

    //Functions
    DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off);
    virtual ~DunePrismAnalyzer(){}; 
    void Analyze();
    void Finalize();
    void SetBranches();
    ////////////////////

    TFile * fin;
    TFile * fout;

    TTree * inTree;
    TTree * gTree;
    TTree * eventTree;
    TTree * gOutTree;
   
    double dimension[3];
    double FV[3];
    double shift;
  
    int nEntries;
    double Enu;
    double EvtVtx[3];
    int nuPID;

    int nstep;
    int ni;
    int PIDi[maxInit];
    double pxi[maxInit];
    double pyi[maxInit];
    double pzi[maxInit];
    double ekini[maxInit];
    double mi[maxInit];

    int PID[maxNQ];
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

    int flagExitBack = 0;
    int flagExitFront = 0;
    int flagExitY = 0;
    int flagExitRight = 0;
    int flagExitLeft = 0;

    double muExitingP[3];

    double eHadOut;
    double eHadIn;
    double eHadTotalDep;

    double eHadTrueCharged;
    double eHadTrueTotal;
    double eMu;

  private:
};

void parse_args(int, char ** );

double detector[3]; 
double fiducialGap[3];
double offset; 
std::string inFileName;
std::string outFileName;



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

