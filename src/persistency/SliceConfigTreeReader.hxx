#ifndef SLICECONFIGTREEREADER_HXX_SEEN
#define SLICECONFIGTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"
#include "TH1D.h"

#include <string>
#include <vector>

struct SliceConfig : public ITreeReader  {

SliceConfig(){}
SliceConfig(std::string const &inputFile, std::string const &inputDir="");

  Double_t XRange[2];
  Double_t Coeff;

  std::string dir;

  std::string TreeName();

  void Reset();
  void Copy(SliceConfig const &);

  void SetBranchAddresses();

  std::vector< std::pair<double,double> > XRanges;
  std::vector<double> Coeffs;

  std::vector<double> XRangeBins;
  std::vector<double> BinCoeffsVector;

  static std::pair< std::vector<double>,std::vector<double> >
    BuildXRangeBinsCoeffs(
      std::vector< std::pair<double, double> > const &XRanges,
      double const *Coeffs, bool SignFlipX=false);

  void ReadTree();

  std::vector< std::pair<double,double> > GetXRanges();
  std::vector<double> GetCoeffs();

  std::vector<double> GetXRangeBins();
  std::vector<double> GetBinCoeffs();

  TH1D *BuildSliceBinningHelper(std::string const &histName);

  static SliceConfig *MakeTreeWriter();

  ~SliceConfig();
};

#endif
