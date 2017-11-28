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

class DunePrismAnalyzer{
  public:

    DunePrismAnalyzer(std::string inFileName, std::string outFileName, double det[3], double fid[3], double off);
    virtual ~DunePrismAnalyzer(){}; 
    void Analyze();
    void Finalize();

    TFile * fin;
    TFile * fout;

    TTree * inTree;
    TTree * gTree;
    TTree * eventTree;
    TTree * gOutTree;
   
    double * detector;
    double * fiducialGap;
    double offset;
  
    int nEntries;
    double Enu;
    double EvtVtx[3];
    int nuPID;
    
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

