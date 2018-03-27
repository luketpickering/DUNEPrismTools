#ifndef FLUXFITRESULTSTREEREADER_HXX_SEEN
#define FLUXFITRESULTSTREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>

struct FluxFitResultsTreeReader {

FluxFitResultsTreeReader();
FluxFitResultsTreeReader(std::string const &treeName, std::string const &inputFile);

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

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(FluxFitResultsTreeReader const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  static FluxFitResultsTreeReader *MakeTreeWriter(TTree *tree, bool IsGaussFit);

  void ReleaseInputFile();

  ~FluxFitResultsTreeReader();
};

#endif
