#ifndef SLICECONFIGTREEREADER_HXX_SEEN
#define SLICECONFIGTREEREADER_HXX_SEEN

#include "TChain.h"
#include "TH1D.h"

#include <string>
#include <vector>

struct SliceConfig {

SliceConfig();
SliceConfig(std::string const &treeName, std::string const &inputFile);

  Double_t XRange[2];
  Double_t Coeff;

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

  std::vector< std::pair<double,double> > XRanges;
  std::vector<double> Coeffs;

  std::vector<double> XRangeBins;
  std::vector<double> BinCoeffsVector;

  static std::pair< std::vector<double>,std::vector<double> >
    BuildXRangeBinsCoeffs(
      std::vector< std::pair<double, double> > const &XRanges,
      double const *Coeffs);

  void ReadTree();

  std::vector< std::pair<double,double> > GetXRanges();
  std::vector<double> GetCoeffs();

  std::vector<double> GetXRangeBins();
  std::vector<double> GetBinCoeffs();

  TH1D *BuildSliceBinningHelper(std::string const &histName);

  static SliceConfig *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~SliceConfig();
};

#endif
