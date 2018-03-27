#ifndef STOPCONFIGTREEREADER_HXX_SEEN
#define STOPCONFIGTREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>

struct StopConfig {

StopConfig();
StopConfig(std::string const &treeName, std::string const &inputFile);

  Double_t ActiveMin[3];
  Double_t ActiveMax[3];
  Double_t VetoGap[3];
  double CenterPosition[3];
  double POTExposure;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(StopConfig const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  static StopConfig *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~StopConfig();
};

#endif
