#ifndef OSCPARAMTREEREADER_HXX_SEEN
#define OSCPARAMTREEREADER_HXX_SEEN

#include "ITreeReader.hxx"

#include <string>

struct OscillationParameters  : public ITreeReader {

OscillationParameters(){}
OscillationParameters(std::string const &inputFile);

  Double_t DipAngle_degrees;
  Double_t OscParams[6];
  Int_t FromNuPDG,ToNuPDG;

  std::string TreeName();

  void Reset();
  void Copy(OscillationParameters const &);

  void SetBranchAddresses();

  static OscillationParameters *MakeTreeWriter();

  ~OscillationParameters();
};

#endif
