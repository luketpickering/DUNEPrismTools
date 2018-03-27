#ifndef SLICECONFIGTREEREADER_HXX_SEEN
#define SLICECONFIGTREEREADER_HXX_SEEN

#include "TChain.h"
#include "TH1D.h"

#include <string>
#include <vector>

struct SliceConfig {

SliceConfig();
SliceConfig(std::string const &treeName, std::string const &inputFile);

  double XRange[2];
  double Coeff;

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(SliceConfig const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  std::vector<double> XRangeBins;
  std::vector<double> CoeffsVector;

  void ReadTree();

  std::vector<double> GetXRangeBins();
  std::vector<double> GetCoeffs();

  TH1D *BuildSliceBinningHelper(std::string const &histName);

  static SliceConfig *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~SliceConfig();
};

#endif
