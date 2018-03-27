#ifndef SIMCONFIGTREEREADER_HXX_SEEN
#define SIMCONFIGTREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>

struct OscillationParameters {

OscillationParameters();
OscillationParameters(std::string const &treeName, std::string const &inputFile);

  Double_t DipAngle_degrees;
  Double_t OscParams[6];

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(OscillationParameters const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  static OscillationParameters *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~OscillationParameters();
};

#endif
