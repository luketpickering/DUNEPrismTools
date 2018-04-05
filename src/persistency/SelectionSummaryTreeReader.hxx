#ifndef SELECTIONSUMMARYTREEREADER_HXX_SEEN
#define SELECTIONSUMMARYTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"

#include <string>

struct SelectionSummary : public ITreeReader  {

SelectionSummary(){}
SelectionSummary(std::string const &inputFile);

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
  Int_t NOOAcceptance;
  Int_t NSel;
  Bool_t SelectOnMuonExit;
  Double_t MuonExitKECut_MeV;
  Double_t HadronicShowerVetoCut_MeV;
  Double_t VertexSelectionFV[3];
  Double_t TotalPOT;

  std::string TreeName();
  void Reset();
  void Copy(SelectionSummary const &);

  void SetBranchAddresses();

  static SelectionSummary *MakeTreeWriter();

  ~SelectionSummary();
};

#endif
