#ifndef FLUXFITRESULTSTREEREADER_HXX_SEEN
#define FLUXFITRESULTSTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"

#include <string>

struct FluxFitResultsTreeReader : public ITreeReader {

FluxFitResultsTreeReader(){}
FluxFitResultsTreeReader(std::string const &inputFile);

  Int_t NFluxes;
  Int_t NIterations;
  Double_t Chi2;
  Double_t RegularisationPenalty;
  Double_t OutOfRangePenalty;
  Double_t FitRange[2];
  Double_t NDOverFDFitScaleFactor;
  Bool_t IsGaussFit;
  Double_t GaussCenter_GeV;
  Double_t GaussWidth_GeV;

  std::string TreeName();

  void Reset();
  void Copy(FluxFitResultsTreeReader const &);

  void SetBranchAddresses();

  static FluxFitResultsTreeReader *MakeTreeWriter(bool IsGaussFit);

  ~FluxFitResultsTreeReader();
};

#endif
