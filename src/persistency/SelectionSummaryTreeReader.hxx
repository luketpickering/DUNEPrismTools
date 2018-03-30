#ifndef SELECTIONSUMMARYTREEREADER_HXX_SEEN
#define SELECTIONSUMMARYTREEREADER_HXX_SEEN

#include "TChain.h"

#include <string>

struct SelectionSummary {

SelectionSummary();
SelectionSummary(std::string const &treeName, std::string const &inputFile);

  Int_t NTotal;
  Int_t NInStops;
  Int_t NNumuCC;
  Int_t NNumuNC;
  Int_t NNumubCC;
  Int_t NNumubNC;
  Int_t NNueCC;
  Int_t NNueNC;
  Int_t NNuebCC;
  Int_t NNuebNC;
  Int_t NOOFV;
  Int_t NSel;
  Bool_t SelectOnMuonExit;
  Double_t MuonExitKECut_MeV;
  Double_t HadronicShowerVetoCut_MeV;
  Double_t VertexSelectionFV[3];

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void Reset();
  void Copy(SelectionSummary const &);

  void SetBranchAddresses();

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  static SelectionSummary *MakeTreeWriter(TTree *tree);

  void ReleaseInputFile();

  ~SelectionSummary();
};

#endif
