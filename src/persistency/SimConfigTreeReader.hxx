#ifndef SIMCONFIGTREEREADER_HXX_SEEN
#define SIMCONFIGTREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>

struct SimConfig {

SimConfig();
SimConfig(std::string const &treeName, std::string const &inputFile);

  Double_t DetMin[3];
  Double_t DetMax[3];
  Double_t FVGap[3];
  Int_t NXSteps;
  Int_t NMaxTrackSteps;
  Double_t POTPerFile;
  Double_t timesep_us;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(SimConfig const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  static SimConfig *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~SimConfig();
};

#endif
