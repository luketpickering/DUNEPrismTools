#ifndef STOPCONFIGTREEREADER_HXX_SEEN
#define STOPCONFIGTREEREADER_HXX_SEEN

#include "BoundingBox.hxx"

#include "ITreeReader.hxx"

#include <string>

struct StopConfig : public ITreeReader   {

StopConfig(){}
StopConfig(std::string const &inputFile);

  Double_t ActiveMin[3];
  Double_t ActiveMax[3];
  Double_t VetoGap[3];
  Double_t CenterPosition[3];
  Double_t POTExposure;

  UInt_t NStops;

  std::string TreeName();

  void Reset();
  void Copy(StopConfig const &);
  void DetermineNStops();

  void SetBranchAddresses();

  std::vector<BoundingBox> GetStopBoundingBoxes(bool RemoveVeto=false,
    std::array<double,3> FVReduction={0,0,0});

  static StopConfig *MakeTreeWriter();

  ~StopConfig();
};

#endif
