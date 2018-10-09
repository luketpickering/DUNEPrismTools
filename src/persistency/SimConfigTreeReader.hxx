#ifndef SIMCONFIGTREEREADER_HXX_SEEN
#define SIMCONFIGTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"

#include <string>

struct SimConfig  : public ITreeReader  {

SimConfig(){}
SimConfig(std::string const &inputFile);

  Double_t DetMin[3];
  Double_t DetMax[3];
  Double_t VetoGap[3];
  Int_t NXSteps;
  Int_t NMaxTrackSteps;
  Double_t POTPerFile;
  Double_t timesep_us;

  std::string TreeName();

  void Reset();
  void Copy(SimConfig const &);

  void SetBranchAddresses();

  static SimConfig *MakeTreeWriter();

  ~SimConfig();
};

#endif
